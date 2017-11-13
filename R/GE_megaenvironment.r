#' Form mega-environments based on winning genotypes from an AMMI-2 model
#'
#' This function fits an AMMI-2 model and then using the fitted values produces
#' a new factor based on the winning genotype in each environment. This new factor is
#' added as a column megaEnv to the input data.
#'
#' @inheritParams GE.AMMI
#'
#' @param method A criterion to determine the best trait, either \code{"max"} or \code{"min"}.
#' @param summaryTable A logical specifying whether a summary table will be returned.
#'
#' @return The input object of class TD with an added extra column megaEnv
#'
#' @examples
#' myDat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env = "env!", genotype = "genotype!", trait = "yld")
#' myTD <- createTD(data = myDat, genotype = "genotype!", env = "env!")
#' TDmegaEnv <- GE.megaEnvironment(TD = myTD, trait = "yld")
#' attr(TDmegaEnv, "summary")
#'
#' @export

GE.megaEnvironment <- function(TD,
                               trait,
                               method = "max",
                               summaryTable = TRUE) {
  #drop factor levels
  TD$genotype <- droplevels(TD$genotype)
  TD$env <- droplevels(TD$env)
  #use AMMI2
  anal <- GE.AMMI(TD = TD, trait = trait, nPC = 2)
  fitted <- anal$fitted
  genoNames <- rownames(fitted)
  envNames  <- colnames(fitted)
  if (method == "max"){
    winPosition <- apply(X = fitted, MARGIN = 2, FUN = which.max)
  } else if (method == "min") {
    winPosition <- apply(X = fitted, MARGIN = 2, FUN = which.min)
  } else {
    stop('Please choose either "max" or "min" as method.')
  }
  winGeno <- genoNames[winPosition]
  winGenoUnique <- unique(winGeno)
  nGenoUnique <- length(winGenoUnique)
  megaLabels <- 1:nGenoUnique
  megaFactor <- factor(winGeno, labels = megaLabels)
  # re-labelling
  levels(megaFactor) <- unique(as.integer(megaFactor))
  # merge factor levels
  megaEnvFactor <- TD$env
  levels(megaEnvFactor) <- as.character(megaFactor)
  TD$megaEnv <- megaEnvFactor
  # re-labelling
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
