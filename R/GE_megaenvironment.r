#' Form mega-environments based on winning genotypes from an AMMI-2 model
#'
#' This function fits an AMMI-2 model and then using the fitted values produces
#' a new factor based on the winning genotype in each environment. This new factor is
#' added as a column megaEnv to the input data. If a column megaEnv already exists the
#' existing data is overwritten with a warning.
#'
#' @inheritParams GeAmmi
#'
#' @param method A criterion to determine the best trait, either \code{"max"} or \code{"min"}.
#' @param summaryTable A logical specifying whether a summary table will be returned.
#'
#' @return The input object of class \code{\link{TD}} with an added extra column megaEnv.
#'
#' @examples
#' data(TDMaize)
#' TDmegaEnv <- GE.megaEnvironment(TD = TDMaize, trait = "yld")
#' attr(TDmegaEnv, "summary")
#'
#' @export

GE.megaEnvironment <- function(TD,
                               trait,
                               method = "max",
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
  if (is.null(method) || !is.character(method) || length(method) > 1 ||
      !method %in% c("min", "max")) {
    stop("method should be either min or max.\n")
  }
  ## drop factor levels
  TD$genotype <- droplevels(TD$genotype)
  TD$env <- droplevels(TD$env)
  ## use AMMI
  AMMI <- GeAmmi(TD = TD, trait = trait, nPC = 2)
  fitted <- AMMI$fitted
  genoNames <- rownames(fitted)
  envNames  <- colnames(fitted)
  if (method == "max") {
    winPosition <- apply(X = fitted, MARGIN = 2, FUN = which.max)
  } else {
    winPosition <- apply(X = fitted, MARGIN = 2, FUN = which.min)
  }
  winGeno <- genoNames[winPosition]
  megaFactor <- factor(winGeno, labels = 1:length(unique(winGeno)))
  ## re-labelling
  levels(megaFactor) <- unique(as.integer(megaFactor))
  ## merge factor levels
  megaEnvFactor <- TD$env
  levels(megaEnvFactor) <- as.character(megaFactor)
  TD$megaEnv <- megaEnvFactor
  ## re-labelling
  levels(TD$megaEnv) <- 1:nlevels(megaEnvFactor)
  if (summaryTable) {
    summTab <- matrix(nrow = length(envNames), ncol = 4)
    colnames(summTab) <- c("Mega env", "env", "genotype", "AMMI estimates")
    megaOrder <- order(megaFactor)
    summTab[, 1] <- megaFactor
    summTab[, 2] <- envNames
    summTab[, 3] <- winGeno
    summTab[, 4] <- sapply(X = 1:length(envNames), FUN = function(x) {
      signif(fitted[winPosition[x], x], 5)
    })
    summTab <- as.data.frame(summTab[megaOrder, ], stringsAsFactors = FALSE)
    summTab[, 4] <- as.numeric(summTab[, 4])
    attr(TD, "summary") <- summTab
  }
  return(TD)
}
