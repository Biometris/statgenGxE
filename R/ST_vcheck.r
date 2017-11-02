#' Identifying outliers
#'
#' Identifies large standardized residuals from a REML analysis.
#'
#' @param data A data frame object containing data used for the mixed modelling analysis.
#' @param stdres A data frame object containing the standardised residuals obtained from the mixed modelling analysis.
#' @param rdf An integer (vector) specifying the residual degrees of freedom.
#' @param rlimit An interger specifying a limit for detection of large standardized residuals;
#' if this is not set, the limit is set automatically according to the number of residual degrees of freedom.
#' @param trait A string (vector) specifying the column name(s) of the trait(s).
#' @param genotype A string specifying the column name of the genotypes.
#' @param entry A string specifying the column name of the entry numbers.
#' @param plotno A string specifying the column name of the plot numbers.
#' @param rep A string specifying the column name of the replicates.
#' @param subblock A string specifying the column name of the sub-blocks.
#' @param row A string specifying the column name of the rows.
#' @param col A string specifying the column name of the columns.
#' @param rowcoordinates A string specifying row coordinates for fitting spatial models. Default, \code{NA}.
#' @param colcoordinates A string specifying col coordinates for fitting spatial models. Default, \code{NA}.
#' @param commonfactor A string vector specifying factors to define similar units; default, \code{commonfactor=genotype}.
#' @param verbose Logical; whether to print a summary table of outliers.
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"), "SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Plot", "Rep", "Row", "Column"),
#'                      trait.names="yield", env ="Env", rowSelect="DRIP05")
#' sb.lmer <- ST.mod.rowcol(Y=mydat, sub.design="res.rowcol", trait="yield",
#'                          genotype="Genotype", rep="Rep", row="Row", col="Column",
#'                          engine="lme4")
#' stdResid0 <- resid(sb.lmer$mfix, scaled=TRUE)
#' rdf0 <- nrow(model.frame(sb.lmer$mfix))-length(fixef(sb.lmer$mfix))
#' vcheck0 <- ST.vcheck(data=mydat, stdres=stdResid0, rdf=rdf0, trait="yield",
#'                      genotype="Genotype", rep="Rep", plotno="Plot",
#'                      row="Row", col="Column", verbose=TRUE, commonfactor=c("Genotype", "Row"))
#'
#' sb.asr <- ST.mod.rowcol(Y=mydat, sub.design="res.rowcol", trait="yield",
#'                         genotype="Genotype", rep="Rep", row="Row", col="Column",
#'                         engine="asreml")
#' stdResid1 <- residuals(sb.asr$mfix, type="stdCond")
#' rdf1 <- sb.asr$mfix$nedf
#' vcheck1 <- ST.vcheck(data=mydat, stdres=stdResid1, rdf=rdf1, trait="yield",
#'                      genotype="Genotype", rep="Rep", plotno="Plot",
#'                      row="Row", col="Column", verbose=TRUE)
#'
#' @export
ST.vcheck <- function(data, stdres, rdf, rlimit=NA, trait, genotype,
                   entry=NA, plotno=NA, rep=NA, subblock=NA, row=NA, col=NA,
                   rowcoordinates=NA, colcoordinates=NA, commonfactor=genotype,
                   verbose = FALSE){
  nr <- nrow(data)
  if(is.vector(stdres)){
    stdres <- data.frame(stdres)
    names(stdres) <- trait
  }
  absresids <- abs(stdres)
  pmat <- data.frame(TraitName=character(0), TraitValue=numeric(0),
                     Genotype=character(0), Entry=character(0), PlotNo=character(0),
                     Replicate=character(0), Subblock=character(0), Row=character(0), Column=character(0),
                     RowPosition=character(0), ColPosition=character(0), Residual=numeric(0), Similar=integer(0))

  testnames <- c(entry, plotno, rep, subblock, row, col, rowcoordinates, colcoordinates)
  na_names <- is.na(testnames)
  testnames <- testnames[!na_names]
  if(!all(testnames %in% names(data))){
    stop(paste(testnames[!(testnames %in% names(data))], collapse=","),
         " not found in the names of data")
  }

  indicator <- data[,trait, drop=F]
  similar <- rep(2, nr)

  if (verbose)
    cat("the residual method are large standardized residuals reported by the mixed model analysis\n")

  if (length(trait)!=length(rdf))
    stop("number of traits does not equal to number of residual degrees of freedom")
  if (ncol(absresids)!=length(rdf))
    stop("number of columns in standardised residuals does not equal to number of residual degrees of freedom")
  if (ncol(absresids)!=length(trait))
    stop("number of traits does not equal to number of columns in standardised residuals")

  for (ii in 1:length(trait)){
    if (is.na(rlimit)){
      if (rdf[ii] <= 20){
        rlimit <- 2
      }else{
        if (rdf[ii] <= 15773){
          rlimit <- qnorm(1 - 0.5/rdf)
        }else{
          rlimit <- 4
        }
      }
    }
    # identifies outliers as large standardized residuals
    tind <- indicator[,trait[ii]] <- (absresids[,trait[ii]] > rlimit)
    noutliers <- sum(tind)
    # list all factors as similar to those with the large residuals
    ncfac <- length(commonfactor)
    geno.rlarge <- ttind <- rep(FALSE, nr)
    for (jj in 1:ncfac){
      ttind <- data[tind,commonfactor[jj]]
      ttind <- droplevels(ttind)
      ttind <- levels(ttind)
      ttind <- data[[commonfactor[jj]]] %in% ttind
      geno.rlarge <- geno.rlarge + ttind
    }
    geno.rlarge <- geno.rlarge>0
    similar[geno.rlarge] <- 1
    similar[tind] <- 0
    nn <- sum(geno.rlarge)
    if (nn>0){
      pmat0 <- data.frame(TraitName=character(nn), TraitValue=numeric(nn),
                          Genotype=character(nn), Entry=character(nn), PlotNo=character(nn),
                          Replicate=character(nn), Subblock=character(nn), Row=character(nn), Column=character(nn),
                          RowPosition=character(nn), ColPosition=character(nn), Residual=numeric(nn), Similar=integer(nn))
      #pmat0$ENVNM <- data[geno.rlarge, env]
      pmat0$TraitValue <- data[geno.rlarge, trait[ii]]
      pmat0$TraitName <- rep(trait[ii], nn)
      pmat0$Genotype <- data[geno.rlarge,genotype]
      pmat0$Residual <- stdres[geno.rlarge,ii]
      pmat0$Similar <- similar[geno.rlarge]
      if (!is.na(entry))
        pmat0$Entry <- data[geno.rlarge,entry]
      else
        pmat0$Entry <- rep(NA,nn)
      if (!is.na(plotno))
        pmat0$PlotNo <- data[geno.rlarge,plotno]
      else
        pmat0$PlotNo <- rep(NA,nn)
      if (!is.na(rep))
        pmat0$Replicate <- data[geno.rlarge,rep]
      else
        pmat0$Replicate <- rep(NA,nn)
      if (!is.na(subblock))
        pmat0$Subblock <- data[geno.rlarge,subblock]
      else
        pmat0$Subblock <- rep(NA,nn)
      if (!is.na(row))
        pmat0$Row <- data[geno.rlarge,row]
      else
        pmat0$Row <- rep(NA,nn)
      if (!is.na(col))
        pmat0$Column <- data[geno.rlarge,col]
      else
        pmat0$Column <- rep(NA,nn)
      if (!is.na(rowcoordinates))
        pmat0$RowPosition <- data[geno.rlarge,rowcoordinates]
      else
        pmat0$RowPosition <- rep(NA,nn)
      if (!is.na(colcoordinates))
        pmat0$ColPosition <- data[geno.rlarge,colcoordinates]
      else
        pmat0$ColPosition <- rep(NA,nn)
      pmat <- rbind(pmat, pmat0)
    }
  }
  if (nrow(pmat)>0 && verbose){
    na_names0 <- c(FALSE, FALSE, FALSE, na_names, FALSE, FALSE)
    print(format(pmat[,!na_names0]), quote = FALSE)
  }
  indicator
}
