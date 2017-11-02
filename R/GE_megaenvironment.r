#' Form mega-environments based on winning genotypes from an AMMI-2 model
#'
#' This function fits an AMMI-2 model and then using the fitted values produces
#' a new factor based on the winning genotype in each environment.
#'
#' @param Y A data frame object.
#' @param trait A character string specifying the name of a trait column
#' @param genotype A character string specifying the name of a genotype column.
#' @param env A character string specifying the name of an environment column.
#' @param megaenv A character string specifying the name of an mega-environment column.
#' This column is to be added into \code{Y}.
#' @param method A criterion to determine the best trait, either \code{"max"} or \code{"min"}.
#' @param summary.table A logical specifying whether a summary table will be returned.
#' @return A data frame object, consisting of Y and a mega-environment factor.
#' @examples
#' mydat <- GE.read.csv(file.path(path.package("RAP"),"F2maize_pheno.csv"),
#'                      env="env!", genotype="genotype!", trait="yld")
#' names(mydat)=c("env", "genotype","yld") 
#' Y <- GE.megaenvironment(Y=mydat, trait="yld", genotype="genotype", 
#'                         env="env", megaenv="megaenv")
#' str(Y)
#' 
#' @export

GE.megaenvironment <- function(Y, trait, genotype, env, megaenv,
method="max", summary.table = T){

  #drop factor levels
  Y[[genotype]] <- droplevels(Y[[genotype]])
  Y[[env]] <- droplevels(Y[[env]])
    
  #use AMMI2
  anal <- GE.AMMI(Y, trait, genotype, env, nPC=2, AMMI2plot=F)
  fitted <- anal$fitted
  
  genonames <- rownames(fitted)
  envnames  <- colnames(fitted)
  
  if (method == "max"){
    win.position <- apply(fitted, 2, which.max)
  }else{
    if (method == "min")
      win.position <- apply(fitted, 2, which.min)
    else
      stop('Please choose either "max" or "min" as method.')
  }
  win.geno <- genonames[win.position]
  win.geno.unique <- unique(win.geno)
  ngeno.unique <- length(win.geno.unique)
  megalabels <- 1:ngeno.unique
  megafactor <-  factor(win.geno, labels=megalabels)
  # re-labelling
  levels(megafactor) <- unique(as.integer(megafactor))
  
  # merge factor levels
  megaEnvFactor <- Y[[env]]
  levels(megaEnvFactor) <- as.character(megafactor)
  Y[[megaenv]] <- megaEnvFactor
  # re-labelling
  levels(Y[[megaenv]]) <- 1:nlevels(megaEnvFactor)
  
  if (summary.table){
    summTab <- matrix(NA, nrow=length(envnames), ncol=4)
    colnames(summTab) <- c("Mega env", env, genotype, "AMMI estimates")
    megaorder <- order(megafactor)
    summTab[,1] <- megafactor
    summTab[,2] <- envnames
    summTab[,3] <- win.geno
    summTab[,4] <- sapply(1:length(envnames), function(x) signif(fitted[win.position[x],x],5))
    summTab <- as.data.frame(summTab[megaorder,], stringsAsFactors = F)
    summTab[,4] <- as.numeric(summTab[,4])
    attr(Y, "summary") <- summTab
  }
  Y
}
