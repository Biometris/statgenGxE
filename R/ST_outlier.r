#' Identifying outliers
#'
#' This function is to find observations that exceed 1.5 times the inerquartile range.
#'
#' @param data A string path where the data list is saved.
#' @param trait A string (vector) specifying the column name(s) of the trait(s).
#' @param genotype A string specifying the column name of the genotypes.
#' @param entry A string specifying the column name of the entries.
#' @param plotNo A string specifying the column name of the plot numbers found.
#' @param rep A string specifying the column name of the replicates.
#' @param subBlock A string specifying the column name of the sub-blocks.
#' @param row A string specifies the column name of the rows.
#' @param col A string specifies the column name of the columns.
#' @param rowCoordinates A string specifying row coordinates for fitting spatial models; default, \code{NA}.
#' @param colCoordinates A string specifying col coordinates for fitting spatial models; default, \code{NA}.
#' @param commonFactor A string vector specifying factors to define similar units; default, \code{commonFactor=genotype}.
#' @param coef this determines how far the plot 'whiskers' extend out from the box.
#' If \code{coef} is positive, the whiskers extend to the most extreme data point which
#' is no more than coef times the length of the box away from the box.
#' A value of zero causes the whiskers to extend to the data extremes (and no outliers be returned).
#' @return a data frame object containing logical values to determine if the observation is an outlier.
#' @examples
#' mydat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames=c("Env","Genotype","Rep","Row"),
#'                      traitNames="yield", env ="Env",rowSelect="HEAT05",
#'                      colSelect=c("Env","Genotype","Rep","yield","Row"))
#' outliers <- ST.outlier(mydat, trait="yield", genotype="Genotype", rep="Rep")
#' #outliers <- ST.outlier(mydat, trait="yield", genotype="Genotype", rep="Rep", coef=1.2)
#' #outliers <- ST.outlier(mydat, trait="yield", genotype="Genotype", row="Row",
#' #                       commonFactor=c("Genotype", "Row"))
#'
#' @import grDevices
#' @export

ST.outlier <- function(data,
                       trait,
                       genotype,
                       entry = NA,
                       plotNo = NA,
                       rep = NA,
                       subBlock = NA,
                       row = NA,
                       col = NA,
                       rowCoordinates = NA,
                       colCoordinates = NA,
                       commonFactor = genotype,
                       coef = 1.5) {
  testNames <- c(trait, genotype, entry, plotNo, rep, subBlock, row, col,
                 rowCoordinates, colCoordinates)
  naNames <- is.na(testNames)
  testNames <- testNames[!naNames]
  if(!all(testNames %in% names(data))){
    stop(paste(testNames[!(testNames %in% names(data))], collapse = ","),
         " not found in the names of data.\n")
  }
  trts <- indicator <- data[, trait, drop = FALSE]
  similar <- rep(2, nrow(data))
  cat(paste("Observations that exceed", coef, "times the inerquartile range\n"))
  pmat <- data.frame(trait = character(0), value = numeric(0), genotype = character(0),
                     entry = character(0), plotNo = character(0), replicate = character(0),
                     subBlock = character(0), row = character(0), column = character(0),
                     rowPosition = character(0), colPosition = character(0), similar = integer(0))
  for (ii in 1:length(trait)) {
    outVals <- boxplot.stats(x = trts[[trait[ii]]], coef = coef)$out
    tInd <- indicator[[trait[ii]]] <- trts[[trait[ii]]] %in% outVals
    if (length(outVals)) {
      ncFac <- length(commonFactor)
      genoOut <- ttInd <- rep(FALSE, nrow(data))
      for (jj in 1:ncFac){
        ttInd <- data[tInd, commonFactor[jj]]
        ttInd <- droplevels(ttInd)
        ttInd <- levels(ttInd)
        ttInd <- data[[commonFactor[jj]]] %in% ttInd
        genoOut <- genoOut + ttInd
      }
      genoOut <- genoOut > 0
      similar[genoOut] <- 1
      similar[tInd] <- 0
      nn <- sum(genoOut)
      pmat0 <- data.frame(trait = character(nn), value = numeric(nn), genotype = character(nn),
                          entry = character(nn), plotNo = character(nn), replicate = character(nn),
                          subBlock = character(nn), row = character(nn), column = character(nn),
                          rowPosition = character(nn), colPosition = character(nn),
                          similar = integer(nn))
      pmat0$value <- data[genoOut, trait[ii]]
      pmat0$trait <- rep(trait[ii], nn)
      pmat0$similar <- similar[genoOut]
      pmat0$genotype <- data[genoOut, genotype]
      if (!is.na(entry)) {
        pmat0$entry <- data[genoOut, entry]
      } else {
        pmat0$entry <- rep(NA, nn)
      }
      if (!is.na(plotNo)) {
        pmat0$plotNo <- data[genoOut, plotNo]
      } else {
        pmat0$plotNo <- rep(NA, nn)
      }
      if (!is.na(rep)) {
        pmat0$replicate <- data[genoOut, rep]
      } else {
        pmat0$replicate <- rep(NA, nn)
      }
      if (!is.na(subBlock)) {
        pmat0$subBlock <- data[genoOut, subBlock]
      } else {
        pmat0$subBlock <- rep(NA, nn)
      }
      if (!is.na(row)) {
        pmat0$row <- data[genoOut, row]
      } else {
        pmat0$row <- rep(NA,nn)
      }
      if (!is.na(col)) {
        pmat0$column <- data[genoOut, col]
      } else {
        pmat0$column <- rep(NA, nn)
      }
      if (!is.na(rowCoordinates)) {
        pmat0$rowPosition <- data[genoOut, rowCoordinates]
      } else {
        pmat0$rowPosition <- rep(NA, nn)
      }
      if (!is.na(colCoordinates)) {
        pmat0$colPosition <- data[genoOut, colCoordinates]
      } else {
        pmat0$colPosition <- rep(NA, nn)
      }
      pmat <- rbind(pmat, pmat0)
    }
  }
  if (nrow(pmat) > 0) {
    naNames0 <- c(FALSE, naNames, FALSE)
    print(format(pmat[, !naNames0]), quote = FALSE)
  }
  invisible(indicator)
}
