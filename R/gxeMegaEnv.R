#' Form mega-environments based on winning genotypes from an AMMI model
#'
#' This function fits an AMMI model and then using the fitted values produces
#' a new factor based on the winning genotype in each environment. This factor
#' is added as a column megaEnv to the input data. If a column megaEnv already
#' exists the existing data is overwritten with a warning.
#'
#' @inheritParams gxeAmmi
#'
#' @param method A character string indicating the criterion to determine
#' the best genotype per environment, either \code{"max"} or \code{"min"}.
#' @param sumTab Should a summary table be added as an attribute to
#' the output and be printed?
#'
#' @return The input object of class \code{\link{TD}} with an added extra
#' column megaEnv.
#'
#' @examples
#' ## Calculate mega-environments for TDMaize and print a summary of the results.
#' TDmegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld")
#' ## Calculate new mega-environments based on the genotypes with the lowest
#' ## value per environment.
#' TDmegaEnv2 <- gxeMegaEnv(TD = TDmegaEnv, trait = "yld", method = "min")
#'
#' @importFrom methods getFunction
#' @export
gxeMegaEnv <- function(TD,
                       trait,
                       method = c("max", "min"),
                       sumTab = TRUE) {
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
  ## Save and then drop factor levels.
  envLevels <- levels(TD$env)
  TD$env <- droplevels(TD$env)
  ## Perform AMMI anlaysis.
  AMMI <- gxeAmmi(TD = TD, trait = trait, nPC = 2)
  fitted <- AMMI$fitted
  ## Extract position of best genotype per environment.
  winPos <- apply(X = fitted, MARGIN = 2,
                  FUN = getFunction(paste0("which.", method)))
  ## Extract best genotype per environment.
  winGeno <- rownames(fitted)[winPos]
  ## Create factor based on best genotypes.
  megaFactor <- factor(winGeno, labels = "")
  ## Merge factor levels to original data.
  TD$megaEnv <- TD$env
  levels(TD$megaEnv) <- as.character(megaFactor)
  ## Reapply saved levels to ensure input and output TD are identical.
  levels(TD$env) <- envLevels
  if (sumTab) {
    ## Create summary table.
    summTab <- data.frame(megaFactor, envNames = colnames(fitted), winGeno,
                          "AMMI estimates" = fitted[matrix(c(winPos, 1:ncol(fitted)),
                                                           ncol = 2)],
                          check.names = FALSE)
    summTab <- summTab[order(megaFactor), ]
    attr(TD, "sumTab") <- summTab
    print(summTab, row.names = FALSE, digits = 5)
  }
  return(TD)
}
