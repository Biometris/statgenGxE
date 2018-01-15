#' Form mega-environments based on winning genotypes from an AMMI model
#'
#' This function fits an AMMI model and then using the fitted values produces
#' a new factor based on the winning genotype in each environment. This factor is
#' added as a column megaEnv to the input data. If a column megaEnv already
#' exists the existing data is overwritten with a warning.
#'
#' @inheritParams gxeAmmi
#'
#' @param method A character string indicating the criterion to determine
#' the best trait, either \code{"max"} or \code{"min"}.
#' @param summaryTable Should a summary table will be printed?
#'
#' @return The input object of class \code{\link{TD}} with an added extra
#' column megaEnv.
#'
#' @examples
#' TDmegaEnv <- gxeMegaEnvironment(TD = TDMaize, trait = "yld")
#' attr(TDmegaEnv, "summary")
#'
#' @export

gxeMegaEnvironment <- function(TD,
                               trait,
                               method = c("max", "min"),
                               summaryTable = TRUE) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (!"env" %in% colnames(TD)) {
    stop("TD should contain a column env to be able to run an AMMI analysis.\n")
  }
  if ("megaEnv" %in% colnames(TD)) {
    warning("TD already contains a column megaEnv. This column will be overwritten.\n")
  }
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TD)) {
    stop("trait has to be a column in TD.\n")
  }
  method <- match.arg(method)
  ## drop factor levels
  TD$genotype <- droplevels(TD$genotype)
  TD$env <- droplevels(TD$env)
  ## Perform AMMI anlaysis.
  AMMI <- gxeAmmi(TD = TD, trait = trait, nPC = 2)
  fitted <- AMMI$fitted
  ## Extract position of best genotype per environment.
  winPosition <- apply(X = fitted, MARGIN = 2,
                       FUN = getFunction(paste0("which.", method)))
  ## Extract best genotype per environment.
  winGeno <- rownames(fitted)[winPosition]
  ## Create factor based on best genotypes.
  megaFactor <- factor(winGeno, labels = "")
  ## Merge factor levels to original data.
  TD$megaEnv <- TD$env
  levels(TD$megaEnv) <- as.character(megaFactor)
  if (summaryTable) {
    ## Create summary table.
    summTab <- data.frame(megaFactor, envNames = colnames(fitted), winGeno,
                          "AMMI estimates" = sapply(X = 1:ncol(fitted),
                                                    FUN = function(x) {
                                                      signif(fitted[winPosition[x], x], 5)
                                                    }),
                          check.names = FALSE)
    summTab <- summTab[order(megaFactor), ]
    attr(TD, "summary") <- summTab
    print(summTab)
  }
  return(TD)
}
