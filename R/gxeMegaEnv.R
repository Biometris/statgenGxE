#' Form mega environments based on fitted values from an AMMI model
#'
#' This function fits an AMMI model and then using the fitted values produces
#' a new factor clustering the trials. This factor is added as a column megaEnv
#' to the input data. If a column megaEnv already exists this column is
#' overwritten with a warning.\cr\cr
#' Mega environments are created by grouping environments based on their best
#' performing genotype; i.e. environments that share the same best genotype
#' belong to the same mega environment.
#'
#' @inheritParams gxeAmmi
#'
#' @param method A character string indicating the criterion to determine
#' the best genotype per environment, either \code{"max"} or \code{"min"}.
#'
#' @return An object of class megaEnv, a list consisting of
#' \describe{
#' \item{TD}{An object of class TD, the TD object used as input to the function
#' with an extra column megaEnv.}
#' \item{summTab}{A data.frame, a summary table containing information on the
#' trials in each mega environment.}
#' \item{trait}{The trait used for calculating the mega environments.}
#' }
#'
#' @examples
#' ## Calculate mega environments for TDMaize.
#' gemegaEnv <- gxeMegaEnv(TD = TDMaize, trait = "yld")
#'
#' ## Calculate new mega environments based on the genotypes with the lowest
#' ## value per environment.
#' gemegaEnv2 <- gxeMegaEnv(TD = TDMaize, trait = "yld", method = "min")
#'
#' @references Atlin, G. N., R. J. Baker, K. B. McRae, and X. Lu. 2000.
#' Selection Response in Subdivided Target Regions. Crop Sci. 40:7-13.
#' \doi{10.2135/cropsci2000.4017}
#'
#' @family mega environments
#'
#' @export
gxeMegaEnv <- function(TD,
                       trials = names(TD),
                       trait,
                       method = c("max", "min"),
                       byYear = FALSE) {
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  trials <- chkTrials(trials, TD)
  TDTot <- Reduce(f = rbind, x = TD[trials])
  chkCol(trait, TDTot)
  chkCol("trial", TDTot)
  chkCol("genotype", TDTot)
  if (byYear) {
    chkCol("year", TDTot)
  }
  method <- match.arg(method)
  if (hasName(x = TDTot, name = "megaEnv")) {
    warning("TD already contains a column megaEnv. This column will",
            "be overwritten.\n", call. = FALSE)
  }
  ## Remove genotypes that contain only NAs
  allNA <- by(TDTot, TDTot[["genotype"]], FUN = function(x) {
    all(is.na(x[trait]))
  })
  TDTot <- TDTot[!TDTot[["genotype"]] %in% names(allNA[allNA]), ]
  rmYear <- FALSE
  if (!byYear) {
    TDTot[["year"]] <- 0
    rmYear <- TRUE
  }
  ## Save and then drop factor levels.
  envLevels <- levels(TDTot[["trial"]])[levels(TDTot[["trial"]]) %in% trials]
  TDTot[["trial"]] <- droplevels(TDTot[["trial"]])
  ## Perform AMMI analysis.
  AMMI <- gxeAmmi(TD = createTD(TDTot), trait = trait, nPC = 2, byYear = byYear)
  ## Extract winning genotype per trial.
  winGeno <- by(data = AMMI$fitted, INDICES = AMMI$fitted[["trial"]],
                FUN = function(trial) {
                  selGeno <- do.call(paste0("which.", method),
                                     args = list(trial[["fittedValue"]]))
                  as.character(trial[["genotype"]])[selGeno]
                })
  ## Extract values for winning genotype per trial.
  winGenoVal <- by(data = AMMI$fitted, INDICES = AMMI$fitted[["trial"]],
                   FUN = function(trial) {
                     do.call(method, args = list(trial[["fittedValue"]]))
                   })
  ## Create factor based on best genotypes.
  megaFactor <- factor(winGeno,
                       labels = paste0("megaEnv_", seq_along(unique(winGeno))))
  ## Merge factor levels to original data.
  TDTot[["megaEnv"]] <- TDTot[["trial"]]
  levels(TDTot[["megaEnv"]]) <- as.character(megaFactor)
  ## Reapply saved levels to ensure input and output TDTot are identical.
  levels(TDTot[["trial"]]) <- envLevels
  ## If year was added, remove if before creating output.
  if (isTRUE(rmYear)) {
    TDTot <- TDTot[-which(colnames(TDTot) == "year")]
  }
  ## Relevel megaEnv so it is in increasing order.
  TDTot[["megaEnv"]] <- factor(as.character(TDTot[["megaEnv"]]))
  TDOut <- createTD(TDTot)
  ## Create summary table.
  summTab <- data.frame("Mega_factor" = megaFactor,
                        Trial = names(winGeno),
                        "Winning_genotype" = as.character(winGeno),
                        "AMMI_estimates" = as.numeric(winGenoVal))
  summTab <- summTab[order(megaFactor), ]
  return(createMegaEnv(TD = TDOut, summTab = summTab, trait = trait))
}
