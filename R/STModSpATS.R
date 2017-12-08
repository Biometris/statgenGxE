#' Fit Single Trial Model using SpATS
#'
#' Fit Single Trial Model using SpATS
#'
#' @inheritParams STRunModel
#'
#' @examples
#' ## Load data.
#' data(TDHeat05)
#'
#' ## Fit basic spatial model - no blocking, no replicates.
#' modSpATS1 <- STModSpATS(TD = TDHeat05, trait = "yield")
#'
#' ## Fit spatial model including replicates - no blocking.
#' modSpATS2 <- STModSpATS(TD = TDHeat05, trait = "yield", design = "res.rowcol")
#'
#' ## Fit spatial model including replicates and blocking.
#' modSpATS3 <- STModSpATS(TD = TDHeat05, trait = "yield", design = "res.ibd")
#'
#' @export

STModSpATS <- function(TD,
                       traits,
                       covariates = NULL,
                       useCheckId = FALSE,
                       design = "rowcol",
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
  if (!is.null(covariates) && (!is.character(covariates) ||
                               !(all(covariates %in% colnames(TD))))) {
    stop("covariates have to be columns in TD.\n")
  }
  for (colName in c("rowCoordinates", "colCoordinates",
                    if (design %in% c("rowcol", "res.rowcol")) c("rowId", "colId"),
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
  ## Compute number of segments.
  nSeg <- c(ceiling(nlevels(TD$colId) / 2), ceiling(nlevels(TD$rowId) / 2))
  ## If valid values for nSeg are provided in control use these instead.
  if ("nSeg" %in% names(control)) {
    nSegCt <- control$nSeg
    if (length(nSegCt) == 1) {
      nSegCt <- rep(x = nSegCt, times = 2)
    }
    if (is.numeric(nSegCt) && length(nSegCt) <= 2 && all(nSegCt >= 1) &&
        all(nSegCt <= c(nlevels(TD$colId), nlevels(TD$rowId)))) {
      nSeg <- nSegCt
    } else {
      warning("Invalid value for control parameter nSeg. Using default values
              instead.\n")
    }
  }
  ## Construct formula for fixed part.
  fixedForm <- as.formula(paste("~",
                                if (useRepIdFix) "repId" else "1",
                                if (useCheckId) "+ checkId",
                                if (!is.null(covariates)) paste(c("", covariates),
                                                                collapse = "+")))
  ## Construct formula for random part. Include repId depending on design.
  if (length(randEff) != 0) {
    randomForm <- as.formula(paste0("~", if (useRepIdFix) "repId:",
                                    "(", paste(randEff, collapse = "+"), ")"))
  } else {
    randomForm <- NULL
  }
  mr <- sapply(X = traits, FUN = function(trait) {
    ## Fit model with genotype random.
    SpATS::SpATS(response = trait, genotype = "genotype",
                 genotype.as.random = TRUE,
                 spatial = ~ SpATS::PSANOVA(colCoordinates, rowCoordinates,
                                            nseg = nSeg, nest.div = c(2, 2)),
                 fixed = fixedForm,
                 random = randomForm,
                 data = TD, control = list(monitoring = 0), ...)
  }, simplify = FALSE)
  mf <- sapply(X = traits, FUN = function(trait) {
    ## Fit model with genotype fixed.
    SpATS::SpATS(response = trait, genotype = "genotype",
                 genotype.as.random = FALSE,
                 spatial = ~ SpATS::PSANOVA(colCoordinates, rowCoordinates,
                                            nseg = nSeg, nest.div = c(2, 2)),
                 fixed = fixedForm,
                 random = randomForm,
                 data = TD, control = list(monitoring = 0), ...)
  }, simplify = FALSE)
  ## Construct SSA object.
  model <- createSSA(mMix = mr, mFix = mf, data = TD, traits = traits,
                     design = design, engine = "SpATS")
  return(model)
}
