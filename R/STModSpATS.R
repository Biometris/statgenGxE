#' Fit Single Trial Model using SpATS
#'
#' Fit Single Trial Model using SpATS
#'
#' @inheritParams ST.run.model
#'
#' @param useCheckId Should a checkId be used as a fixed parameter in the model?\cr
#' If \code{TRUE} \code{TD} has to contain a column 'checkId'.
#'
#' @examples
#' ## Load data
#' data(TDHeat05)
#'
#' ## Fit basic spatial model - no blocking, no replicates
#' modSpATS1 <- STModSpATS(TD = TDHeat05, trait = "yield")
#'
#' ## Fit spatial model including replicates - no blocking
#' modSpATS2 <- STModSpATS(TD = TDHeat05, trait = "yield", design = "res.rowcol")
#'
#' ## Fit spatial model including replicates and blocking
#' modSpATS3 <- STModSpATS(TD = TDHeat05, trait = "yield", design = "res.ibd")
#'
#' @export

STModSpATS <- function(TD,
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
  for (colName in c("rowId", "colId", "rowCoordinates", "colCoordinates",
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
  useRepIdFix <- design %in% c("res.ibd", "res.rowcol")
  ## Indicate extra random effects.
  if (design %in% c("ibd", "res.ibd")) {
    randEff <- "subBlock"
  } else if (design %in% c("rowcol", "res.rowcol")) {
    randEff <- c("rowId", "colId")
  } else if (design == "rcbd") {
    randEff <- "repId"
  }
  ## Compute number of segments.
  nSeg <- c(ceiling(nlevels(TD$colId) / 2), ceiling(nlevels(TD$rowId) / 2))
  ## Construct formula for fixed part.
  fixedForm <- as.formula(paste("~",
                                if (useRepIdFix) "repId" else "1",
                                if (useCheckId) "+ checkId",
                                if (!is.null(covariates)) paste(c("", covariates),
                                                               collapse = "+")))
  ## Construct formula for random part. Include repId depending on design.
  randomForm <- as.formula(paste0("~", if (useRepIdFix) "repId:",
                                  "(", paste(randEff, collapse = "+"), ")"))
  ## Fit model with genotype random.
  mr <- SpATS::SpATS(response = trait, genotype = "genotype",
                     genotype.as.random = TRUE,
                     spatial = ~ SpATS::PSANOVA(colCoordinates, rowCoordinates, nseg = nSeg),
                     fixed = fixedForm,
                     random = randomForm,
                     data = TD, control = list(monitoring = 2), ...)
  ## Fit model with genotype fixed.
  mf <- SpATS::SpATS(response = trait, genotype = "genotype",
                     genotype.as.random = FALSE,
                     spatial = ~ SpATS::PSANOVA(colCoordinates, rowCoordinates, nseg = nSeg),
                     fixed = fixedForm,
                     random = randomForm,
                     data = TD, control = list(monitoring = 2), ...)
  ## Construct SSA object.
  model = createSSA(mMix = mr, mFix = mf, data = TD, trait = trait,
                    genotype = "genotype",
                    repId = if (design == "res.ibd") "repId" else NULL,
                    design = design, engine = "SpATS")
  return(model)
}
