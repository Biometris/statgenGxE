#' Form mega-environments based on winning genotypes from an AMMI-2 model
#'
#' This function fits an AMMI-2 model and then using the fitted values produces
#' a new factor based on the winning genotype in each environment.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying the name of a trait column
#' @param genotype A character string specifying the name of a genotype column.
#' @param env A character string specifying the name of an environment column.
#' @param megaEnv A character string specifying the name of an mega-environment column.
#' This column is to be added into \code{Y}.
#' @param method A criterion to determine the best trait, either \code{"max"} or \code{"min"}.
#' @param summaryTable A logical specifying whether a summary table will be returned.
#' @return A data frame object, consisting of Y and a mega-environment factor.
#' @examples
#' mydat <- GE.read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld")
#' Y <- GE.megaEnvironment(Y=mydat, trait="yld", genotype="genotype",
#'                         env="env", megaEnv="megaEnv")
#' str(Y)
#'
#' @export

GE.megaEnvironment <- function(Y,
                               trait,
                               genotype,
                               env,
                               megaEnv,
                               method = "max",
                               summaryTable = TRUE) {
  #drop factor levels
  Y[[genotype]] <- droplevels(Y[[genotype]])
  Y[[env]] <- droplevels(Y[[env]])
  #use AMMI2
  anal <- GE.AMMI(Y, trait, genotype, env, nPC = 2, AMMI2plot = FALSE)
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
  megaEnvFactor <- Y[[env]]
  levels(megaEnvFactor) <- as.character(megaFactor)
  Y[[megaEnv]] <- megaEnvFactor
  # re-labelling
  levels(Y[[megaEnv]]) <- 1:nlevels(megaEnvFactor)
  if (summaryTable){
    summTab <- matrix(nrow = length(envNames), ncol = 4)
    colnames(summTab) <- c("Mega env", env, genotype, "AMMI estimates")
    megaOrder <- order(megaFactor)
    summTab[, 1] <- megaFactor
    summTab[, 2] <- envNames
    summTab[, 3] <- winGeno
    summTab[, 4] <- sapply(X = 1:length(envNames), FUN = function(x) {
      signif(fitted[winPosition[x], x], 5)
    })
    summTab <- as.data.frame(summTab[megaOrder, ], stringsAsFactors = FALSE)
    summTab[, 4] <- as.numeric(summTab[, 4])
    attr(Y, "summary") <- summTab
  }
  return(Y)
}
