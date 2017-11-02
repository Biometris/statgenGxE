#' Identifying outliers
#'
#' This function is to find observations that exceed 1.5 times the inerquartile range.
#'
#' @param data A string path where the data list is saved.
#' @param trait A string (vector) specifying the column name(s) of the trait(s).
#' @param genotype A string specifying the column name of the genotypes.
#' @param entry A string specifying the column name of the entries.
#' @param plotno A string specifying the column name of the plot numbers found.
#' @param rep A string specifying the column name of the replicates.
#' @param subblock A string specifying the column name of the sub-blocks.
#' @param row A string specifies the column name of the rows.
#' @param col A string specifies the column name of the columns.
#' @param rowcoordinates A string specifying row coordinates for fitting spatial models; default, \code{NA}.
#' @param colcoordinates A string specifying col coordinates for fitting spatial models; default, \code{NA}.
#' @param commonfactor A string vector specifying factors to define similar units; default, \code{commonfactor=genotype}.
#' @param coef this determines how far the plot 'whiskers' extend out from the box.
#' If \code{coef} is positive, the whiskers extend to the most extreme data point which
#' is no more than coef times the length of the box away from the box.
#' A value of zero causes the whiskers to extend to the data extremes (and no outliers be returned).
#' @return a data frame object containing logical values to determine if the observation is an outlier.
#' @examples
#' mydat <- ST.read.csv(file.path(path.package("RAP"),"SB_yield.csv"),
#'                      factor.names=c("Env","Genotype","Rep","Row"),
#'                      trait.names="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","yield","Row"))
#' outliers <- ST.outlier(mydat, trait="yield", genotype="Genotype", rep="Rep")
#' #outliers <- ST.outlier(mydat, trait="yield", genotype="Genotype", rep="Rep", coef=1.2)
#' #outliers <- ST.outlier(mydat, trait="yield", genotype="Genotype", row="Row",
#' #                       commonfactor=c("Genotype", "Row"))
#' 
#' @export


ST.outlier <- function(data, trait, genotype, entry=NA, plotno=NA, rep=NA, subblock=NA, row=NA, col=NA,
              rowcoordinates=NA, colcoordinates=NA, commonfactor=genotype, coef = 1.5){

  testnames <- c(trait, genotype, entry, plotno, rep, subblock, row, col, rowcoordinates, colcoordinates)
  na_names <- is.na(testnames)
  testnames <- testnames[!na_names]
  if(!all(testnames %in% names(data))){
    stop(paste(testnames[!(testnames %in% names(data))],collapse=","),
    " not found in the names of data")
  }
  trts <- indicator <- data[,trait, drop=F]
  similar <- rep(2, nrow(data))
  cat(paste("Observations that exceed", coef,"times the inerquartile range\n"))
  pmat <- data.frame(Trait=character(0), Value=numeric(0), Genotype=character(0), Entry=character(0), PlotNo=character(0), Replicate=character(0),
  Subblock=character(0), Row=character(0), Column=character(0), RowPosition=character(0), ColPosition=character(0), Similar=integer(0))
  for (ii in 1:length(trait)){
    outvals <- boxplot.stats(trts[[trait[ii]]], coef=coef)$out
    tind <- indicator[[trait[ii]]] <- trts[[trait[ii]]]%in%outvals
    if (length(outvals)){
      ncfac <- length(commonfactor)
      geno.out <- ttind <- rep(FALSE, nrow(data))
      for (jj in 1:ncfac){
        ttind <- data[tind,commonfactor[jj]]
        ttind <- droplevels(ttind)
        ttind <- levels(ttind)
        ttind <- data[[commonfactor[jj]]] %in% ttind
        geno.out <- geno.out + ttind
      }
      geno.out <- geno.out>0
      similar[geno.out] <- 1
      similar[tind] <- 0
      nn <- sum(geno.out)
      pmat0 <- data.frame(Trait=character(nn), Value=numeric(nn), Genotype=character(nn), Entry=character(nn), PlotNo=character(nn), Replicate=character(nn),
      Subblock=character(nn), Row=character(nn), Column=character(nn), RowPosition=character(nn), ColPosition=character(nn), Similar=integer(nn))
      pmat0$Value <- data[geno.out, trait[ii]]
      pmat0$Trait <- rep(trait[ii], nn)
      pmat0$Similar <- similar[geno.out]
      pmat0$Genotype <- data[geno.out,genotype]
      if (!is.na(entry))
        pmat0$Entry <- data[geno.out,entry]
      else
        pmat0$Entry <- rep(NA,nn)
      if (!is.na(plotno))
        pmat0$PlotNo <- data[geno.out,plotno]
      else
        pmat0$PlotNo <- rep(NA,nn)
      if (!is.na(rep))
        pmat0$Replicate <- data[geno.out,rep]
      else
        pmat0$Replicate <- rep(NA,nn)
      if (!is.na(subblock))
        pmat0$Subblock <- data[geno.out,subblock]
      else
        pmat0$Subblock <- rep(NA,nn)
      if (!is.na(row))
        pmat0$Row <- data[geno.out,row]
      else
        pmat0$Row <- rep(NA,nn)
      if (!is.na(col))
        pmat0$Column <- data[geno.out,col]
      else
        pmat0$Column <- rep(NA,nn)
      if (!is.na(rowcoordinates))
        pmat0$RowPosition <- data[geno.out,rowcoordinates]
      else
        pmat0$RowPosition <- rep(NA,nn)
      if (!is.na(colcoordinates))
        pmat0$ColPosition <- data[geno.out,colcoordinates]
      else
        pmat0$ColPosition <- rep(NA,nn)
      pmat <- rbind(pmat, pmat0)
    }
  }
  if (nrow(pmat)>0){
    na_names0 <- c(FALSE, na_names, FALSE)
    print(format(pmat[,!na_names0]), quote = FALSE)
  }
  invisible(indicator)
}
