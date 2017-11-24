#' Identifying outliers
#'
#' This function is to find observations that exceed \code{coef} times the interquartile range.
#'
#' @inheritParams ST.vcheck
#'
#' @param coef this determines how far the plot 'whiskers' extend out from the box.
#' If \code{coef} is positive, the whiskers extend to the most extreme data point which
#' is no more than coef times the length of the box away from the box.
#' A value of zero causes the whiskers to extend to the data extremes
#' (and no outliers be returned).
#'
#' @return a data.frame object containing logical values to determine if the observation is
#' an outlier.
#'
#' @examples
#' myDat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames = c("Env", "Genotype", "Rep", "Row"),
#'                      traitNames = "yield", env = "Env", rowSelect = "HEAT05",
#'                      colSelect = c("Env", "Genotype", "Rep", "yield", "Row"))
#' myTD <- createTD(data = myDat, genotype = "Genotype", env = "Env")
#' outliers <- ST.outlier(TD = myTD, traits = "yield", repId = "Rep")
#'
#' @import grDevices
#' @export

ST.outlier <- function(TD,
                       traits,
                       entry = NA,
                       plotNo = NA,
                       repId = NULL,
                       subBlock = NULL,
                       rowId = NULL,
                       colId = NULL,
                       rowCoordinates = NULL,
                       colCoordinates = NULL,
                       commonFactor = "genotype",
                       coef = 1.5) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (is.null(traits) || !is.character(traits) || !all(traits %in% colnames(TD))) {
    stop("trait has to be a vector of columns in TD.\n")
  }
  for (param in c(repId, subBlock, rowId, colId, rowCoordinates,
                  colCoordinates, commonFactor)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !param %in% colnames(TD))) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  if (is.null(coef) || !is.numeric(coef) || length(coef) > 1 || coef < 0) {
    stop("coef should be a positive numerical value.\n")
  }
  trts <- indicator <- TD[, traits, drop = FALSE]
  similar <- rep(2, nrow(TD))
  cat(paste("Observations that exceed", coef, "times the interquartile range\n"))
  pMat <- data.frame(trait = character(0), value = numeric(0), genotype = character(0),
                     entry = character(0), plotNo = character(0), replicate = character(0),
                     subBlock = character(0), rowId = character(0), column = character(0),
                     rowPosition = character(0), colPosition = character(0), similar = integer(0))
  for (ii in 1:length(traits)) {
    outVals <- boxplot.stats(x = trts[[traits[ii]]], coef = coef)$out
    tInd <- indicator[[traits[ii]]] <- trts[[traits[ii]]] %in% outVals
    if (length(outVals)) {
      ncFac <- length(commonFactor)
      genoOut <- ttInd <- rep(FALSE, nrow(TD))
      for (jj in 1:ncFac) {
        ttInd <- TD[tInd, commonFactor[jj]]
        ttInd <- droplevels(ttInd)
        ttInd <- levels(ttInd)
        ttInd <- TD[[commonFactor[jj]]] %in% ttInd
        genoOut <- genoOut + ttInd
      }
      genoOut <- genoOut > 0
      similar[genoOut] <- 1
      similar[tInd] <- 0
      nn <- sum(genoOut)
      pMat0 <- data.frame(trait = character(nn), value = numeric(nn), genotype = character(nn),
                          entry = character(nn), plotNo = character(nn), replicate = character(nn),
                          subBlock = character(nn), rowId = character(nn), column = character(nn),
                          rowPosition = character(nn), colPosition = character(nn),
                          similar = integer(nn))
      pMat0$value <- TD[genoOut, traits[ii]]
      pMat0$trait <- rep(traits[ii], nn)
      pMat0$similar <- similar[genoOut]
      pMat0$genotype <- TD[genoOut, "genotype"]
      if (!is.na(entry)) {
        pMat0$entry <- TD[genoOut, entry]
      } else {
        pMat0$entry <- rep(NA, nn)
      }
      if (!is.na(plotNo)) {
        pMat0$plotNo <- TD[genoOut, plotNo]
      } else {
        pMat0$plotNo <- rep(NA, nn)
      }
      if (!is.null(repId)) {
        pMat0$replicate <- TD[genoOut, repId]
      } else {
        pMat0$replicate <- rep(NA, nn)
      }
      if (!is.null(subBlock)) {
        pMat0$subBlock <- TD[genoOut, subBlock]
      } else {
        pMat0$subBlock <- rep(NA, nn)
      }
      if (!is.null(rowId)) {
        pMat0$rowId <- TD[genoOut, rowId]
      } else {
        pMat0$rowId <- rep(NA, nn)
      }
      if (!is.null(colId)) {
        pMat0$column <- TD[genoOut, colId]
      } else {
        pMat0$column <- rep(NA, nn)
      }
      if (!is.null(rowCoordinates)) {
        pMat0$rowPosition <- TD[genoOut, rowCoordinates]
      } else {
        pMat0$rowPosition <- rep(NA, nn)
      }
      if (!is.null(colCoordinates)) {
        pMat0$colPosition <- TD[genoOut, colCoordinates]
      } else {
        pMat0$colPosition <- rep(NA, nn)
      }
      pMat <- rbind(pMat, pMat0)
    }
  }
  if (nrow(pMat) > 0) {
  #  naNames0 <- c(FALSE, naNames, FALSE)
   # print(format(pMat[, !naNames0]), quote = FALSE)
  }
  invisible(indicator)
}
