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
                       trials = names(TD),
                       trait,
                       method = c("max", "min"),
                       sumTab = TRUE) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (!is.character(trials) || !all(trials %in% names(TD))) {
    stop("All trials should be in TD.")
  }
  TDTot <- Reduce(f = rbind, x = TD[trials])
  if (is.null(trait) || !is.character(trait) || length(trait) > 1 ||
      !trait %in% colnames(TDTot)) {
    stop("trait has to be a column in TD.\n")
  }
  if (!"trial" %in% colnames(TDTot)) {
    stop("TD should contain a column trial to be able to run an AMMI analysis.\n")
  }
  if ("megaEnv" %in% colnames(TDTot)) {
    warning("TD already contains a column megaEnv. This column will be overwritten.\n")
  }
  method <- match.arg(method)
  ## Save and then drop factor levels.
  envLevels <- levels(TDTot$trial)
  TDTot$trial <- droplevels(TDTot$trial)
  ## Perform AMMI anlaysis.
  AMMI <- gxeAmmi(TD = TD, trait = trait, nPC = 2)
  fitted <- AMMI$fitted
  ## Extract position of best genotype per trial.
  winPos <- apply(X = fitted, MARGIN = 2,
                  FUN = getFunction(paste0("which.", method)))
  ## Extract best genotype per trial.
  winGeno <- rownames(fitted)[winPos]
  ## Create factor based on best genotypes.
  megaFactor <- factor(winGeno, labels = "")
  ## Merge factor levels to original data.
  TDTot$megaEnv <- TDTot$trial
  levels(TDTot$megaEnv) <- as.character(megaFactor)
  ## Reapply saved levels to ensure input and output TDTot are identical.
  levels(TDTot$trial) <- envLevels
  TDTot <- createTD(TDTot)
  if (sumTab) {
    ## Create summary table.
    summTab <- data.frame(megaFactor, envNames = colnames(fitted), winGeno,
                          "AMMI estimates" = fitted[matrix(c(winPos, 1:ncol(fitted)),
                                                           ncol = 2)],
                          check.names = FALSE)
    summTab <- summTab[order(megaFactor), ]
    attr(TDTot, "sumTab") <- summTab
    print(summTab, row.names = FALSE, digits = 5)
  }
  return(TDTot)
}
