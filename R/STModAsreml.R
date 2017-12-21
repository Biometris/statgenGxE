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
  if (!is.character(trySpatial)) {
    #if (is.character(trySpatial)) {
    #  stop("Spatial models not yet implemented for SpATS.\n")
    #}
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
    ## Create tempfile to suppress asreml output messages.
    tmp <- tempfile()
    ## Construct formula for fixed part.
    fixedForm <- paste("~",
                       if (useRepIdFix) "repId" else "1",
                       if (useCheckId) "+ checkId",
                       if (!is.null(covariates)) paste(c("", covariates),
                                                       collapse = "+"))
    ## Construct formula for random part. Include repId depending on design.
    if (length(randEff) != 0) {
      randomfTraitorm <- paste0(if (useRepIdFix) "repId:",
                                paste(randEff, collapse = paste("+", if (useRepIdFix) "repId:")))
    } else {
      randomfTraitorm <- character()
    }
    ## Create empty base lists.
    mr <- mfTrait <- setNames(vector(mode = "list", length = length(traits)),
                              traits)
    for (trait in traits) {
      if ("random" %in% what) {
        ## Fit model with genotype random.
        sink(file = tmp)
        mrTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm)),
                                  random = as.formula(paste("~", randomfTraitorm,
                                                            if (length(randomfTraitorm) != 0) "+",
                                                            "genotype")),
                                  rcov = ~ units, aom = TRUE, data = TD, ...)
        sink()
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
        ## evaluate call terms in mr and mfTrait so predict can be run.
        mrTrait$call$fixed <- eval(mrTrait$call$fixed)
        mrTrait$call$random <- eval(mrTrait$call$random)
        mrTrait$call$rcov <- eval(mrTrait$call$rcov)
        # Run predict.
        mrTrait <- predictAsreml(mrTrait, TD = TD)
        mrTrait$call$data <- substitute(TD)
        mr[[trait]] <- mrTrait
      }
      if ("fixed" %in% what) {
        ## Fit model with genotype fixed.
        if (!"random" %in% what) {
          GParamTmp <- NULL
        }
        sink(file = tmp)
        if (length(randomfTraitorm) != 0) {
          mfTraitTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm,
                                                                  "+ genotype")),
                                         random = as.formula(paste("~", randomfTraitorm)),
                                         rcov = ~ units, G.param = GParamTmp, aom = TRUE,
                                         data = TD, ...)
        } else {
          mfTraitTrait <- asreml::asreml(fixed = as.formula(paste(trait, fixedForm,
                                                                  "+ genotype")),
                                         rcov = ~ units, G.param = GParamTmp, aom = TRUE,
                                         data = TD, ...)
        }
        sink()
        mfTraitTrait$call$fixed <- eval(mfTraitTrait$call$fixed)
        mfTraitTrait$call$random <- eval(mfTraitTrait$call$random)
        mfTraitTrait$call$rcov <- eval(mfTraitTrait$call$rcov)
        ## Construct assocForm for use in associate in predict.
        if (useCheckId) {
          assocForm <- as.formula("~ checkId:genotype")
        } else {
          assocForm <- as.formula("~ NULL")
        }
        ## Run predict.
        mfTraitTrait <- predictAsreml(mfTraitTrait, TD = TD, associate = assocForm)
        mfTraitTrait$call$data <- substitute(TD)
        mfTrait[[trait]] <- mfTraitTrait
      }
    }
    unlink(tmp)
    ## Construct SSA object.
    model <- createSSA(mRand = if ("random" %in% what) mr else NULL,
                       mFix = if ("fixed" %in% what) mfTrait else NULL,
                       data = TD, traits = traits,
                       design = design, engine = "asreml")
  } else {
    model <- bestSpatMod(TD = TD, traits = traits, trySpatial = trySpatial,
                         criterion = "AIC", useCheckId = useCheckId,
                         design = design, covariates = covariates, ...)
  }
  return(model)
}

#' Helper function for calculating best spatial model using asreml.
#' @keywords internal
bestSpatMod <- function(TD,
                        traits,
                        trySpatial = "always",
                        criterion = "AIC",
                        useCheckId = FALSE,
                        design = "rowcol",
                        covariates = NULL,
                        ...) {
  ## Create tempfile to suppress asreml output messages.
  tmp <- tempfile()
  useRepIdFix <- design == "res.rowcol"
  ## Define random terms of models to try.
  randomTerm <- c(rep(x = c("NULL", "units"), each = 3),
                  "repId:rowId", "repId:colId", "repId:rowId + repId:colId",
                  "repId:rowId + units",
                  "repId:colId + units",
                  "repId:rowId + repId:colId + units")
  if (!useRepIdFix) {
    ## If no repId remove this from randomTerm
    randomTerm <- gsub(pattern = "repId:", replacement = "", x = randomTerm)
  }
  if (trySpatial == "ifregular") {
    ## Define spatial terms of models to try.
    spatialChoice <- rep(x = c("AR1(x)id", "id(x)AR1", "AR1(x)AR1"), times = 4)
    spatialTerm <- rep(x = c("ar1(rowId):colId",
                             "rowId:ar1(colId)",
                             "ar1(rowId:ar1(colId)"),
                       times = 4)
  } else if (trySpatial == "always") {
    spatialChoice <- rep(x = c("exp(x)id", "id(x)exp",
                               "isotropic exponential"), times = 4)
    spatialTerm <- rep(x = c("exp(rowCoordinates):colCoordinates",
                             "rowCoordinates:exp(colCoordinates)",
                             "iexp(rowCoordinates,colCoordinates)"),
                       times = 4)
  }
  ## Create empty base lists.
  mr <- mf <- spatial <- setNames(vector(mode = "list", length = length(traits)),
                                  traits)
  for (trait in traits) {
    ## Create formula for the fixed part.
    fixedFormR <- as.formula(paste(trait, "~",
                                   if (useRepIdFix) "repId" else "1",
                                   if (useCheckId) "+ checkId",
                                   if (!is.null(covariates)) paste(c("", covariates),
                                                                   collapse = "+")))
    ## Fit model with genotype random for all different random/spatial terms.
    for (i in 1:length(randomTerm)) {
      sink(file = tmp)
      mrTrait <- asreml::asreml(fixed = fixedFormR,
                                random = as.formula(paste("~ genotype +",
                                                          randomTerm[i])),
                                rcov = as.formula(paste("~", spatialTerm[i])),
                                aom = TRUE, data = TD, ...)
      sink()
      ## If current model is better than best so far based on chosen criterion
      ## define best model as current model.
      if (i == 1) {
        bestModelTrait <- mrTrait
        bestLoc <- 1
      } else {
        if (criterion == "AIC") {
          criterionCur  <- -2 * mrTrait$loglik + 2 * length(mrTrait$gammas)
          criterionPrev <- -2 * bestModelTrait$loglik +
            2 * length(bestModelTrait$gammas)
        } else {
          criterionCur  <- -2 * mrTrait$loglik +
            log(length(mrTrait$fitted.values)) * length(mrTrait$gammas)
          criterionPrev <- -2 * bestModelTrait$loglik +
            log(length(bestModelTrait$fitted.values)) * length(bestModelTrait$gammas)
        }
        if (criterionCur < criterionPrev) {
          bestModelTrait <- mrTrait
          bestLoc <- i
        }
      }
    }
    fixedFormfTrait <- as.formula(paste(deparse(fixedFormR), "+ genotype"))
    ## Constrain variance of the variance components to be fixed as the values in the best model.
    GParamTmp <- mrTrait$G.param
    for (randEf in c("rowId", "colId")) {
      ## When there are no replicates the structure is [[randEf]][[randEf]]
      ## otherwise it is [[repId:randEf]][[repId]]
      GParamTmp[[paste0(ifelse(useRepIdFix, "repId:", ""),
                        randEf)]][[ifelse(useRepIdFix, "repId", randEf)]]$con <- "F"
    }
    sink(file = tmp)
    ## Fit the model with genotype fixed only for the best model.
    mfTrait <- asreml::asreml(fixed = fixedFormfTrait,
                              random = as.formula(paste("~", randomTerm[bestLoc])),
                              rcov = as.formula(paste("~", spatialTerm[bestLoc])),
                              G.param = GParamTmp, aom = TRUE, data = TD, ...)
    sink()
    ## evaluate call terms in bestNidek and mfTrait so predict can be run.
    bestModelTrait$call$fixed <- eval(bestModelTrait$call$fixed)
    bestModelTrait$call$random <- eval(bestModelTrait$call$random)
    bestModelTrait$call$rcov <- eval(bestModelTrait$call$rcov)
    # Run predict.
    bestModelTrait <- predictAsreml(bestModelTrait, TD = TD)
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
    mfTrait <- predictAsreml(mfTrait, TD = TD, associate = assocForm)
    mr[[trait]] <- bestModelTrait
    mf[[trait]] <- mfTrait
    spatial[[trait]] <- spatialChoice[bestLoc]
  }
  unlink(tmp)
  model <- createSSA(mRand = mr,
                     mFix = mf,
                     data = TD,
                     traits = trait,
                     design = design,
                     spatial = spatial,
                     engine = "asreml")
  return(model)
}




