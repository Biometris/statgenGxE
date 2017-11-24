#' Identifying outliers
#'
#' Identifies large standardized residuals from a REML analysis.
#'
#' @inheritParams ST.run.model
#'
#' @param traits A character vector specifying the selected traits.
#' @param stdRes A vector or data.frame object containing the standardised residuals obtained
#' from the mixed modelling analysis.
#' @param rDf An integer (vector) specifying the residual degrees of freedom.
#' @param rLimit An integer specifying a limit for detection of large standardized residuals;
#' if this is not set, the limit is set automatically according to the number of residual
#' degrees of freedom.
#' @param entry A string specifying the column name of the entry numbers.
#' @param plotNo A string specifying the column name of the plot numbers.
#' @param commonFactor A string vector specifying factors to define similar units; default,
#' \code{commonFactor = "genotype"}.
#' @param verbose Logical; whether to print a summary table of outliers.
#'
#' @examples
#' myDat <- ST.read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
#'                      factorNames = c("Env", "Genotype", "Rep", "Plot", "Row", "Column"),
#'                      traitNames = "yield", env = "Env", rowSelect = "HEAT05",
#'                      colSelect = c("Env", "Genotype", "Rep", "Row", "Column",
#'                                    "Plot", "yield"))
#' myTD <- createTD(data = myDat, genotype = "Genotype", env = "Env")
#' myModel <- ST.mod.rowcol(TD = myTD, subDesign = "res.rowcol", trait = "yield",
#'                          repId = "Rep", rowId = "Row", colId = "Column",
#'                          engine = "lme4") #engine = "asreml"
#' stdResid <- resid(myModel$mFix, scaled = TRUE)
#' rDf <- nrow(model.frame(myModel$mFix)) - length(lme4::fixef(myModel$mFix))
#' vCheck <- ST.vcheck(TD = myTD, stdRes = stdResid, rDf = rDf, traits = "yield",
#'                     repId = "Rep", plotNo = "Plot", rowId = "Row", colId = "Column",
#'                     verbose = TRUE, commonFactor = c("genotype", "Row"))
#'
#' @export
ST.vcheck <- function(TD,
                      stdRes,
                      rDf,
                      rLimit = NULL,
                      traits,
                      entry = NULL,
                      plotNo = NULL,
                      repId = NULL,
                      subBlock = NULL,
                      rowId = NULL,
                      colId = NULL,
                      rowCoordinates = NULL,
                      colCoordinates = NULL,
                      commonFactor = "genotype",
                      verbose = FALSE) {
  ## Checks.
  if (missing(TD) || !inherits(TD, "TD")) {
    stop("TD should be a valid object of class TD.\n")
  }
  if (missing(stdRes) || (!is.data.frame(stdRes) && !is.vector(stdRes)) ||
      !is.numeric(stdRes)) {
    stop("stdRes should be a data.frame or vector containing numerical values.\n")
  }
  if (missing(rDf) || !is.vector(rDf) || !is.numeric(rDf) ||
      any(rDf != round(rDf)) || any(rDf < 1)) {
    stop("rDf should be a vector containing positive integers.\n")
  }
  if (!is.null(rLimit) && (!is.numeric(rLimit) || length(rLimit) > 1 ||
      rLimit < 0 || rLimit != round(rLimit))) {
    stop("rLimit should be NULL or a positive numerical value.\n")
  }
  if (is.null(traits) || !is.character(traits) || !all(traits %in% colnames(TD))) {
    stop("trait has to be a vector of columns in TD.\n")
  }
  for (param in c(entry, plotNo, repId, subBlock, rowId, colId, rowCoordinates,
                  colCoordinates, commonFactor)) {
    if (!is.null(param) && (!is.character(param) || length(param) > 1 ||
                            !param %in% colnames(TD))) {
      stop(paste(deparse(param), "has to be NULL or a column in data.\n"))
    }
  }
  if (length(traits) != length(rDf)) {
    stop("number of traits should be equal to the number of residual degrees of freedom.\n")
  }
  if (is.vector(stdRes)) {
    stdRes <- data.frame(stdRes)
    names(stdRes) <- traits
  }
  absResids <- abs(stdRes)
  if (ncol(absResids) != length(rDf)) {
    stop("number of columns in standardised residuals should be equal to number of residual
         degrees of freedom.\n")
  }
  pMat <- data.frame(traitName = character(0), traitValue = numeric(0),
                     genotype = character(0), entry = character(0), plotNo = character(0),
                     replicate = character(0), subBlock = character(0), rowId = character(0),
                     column = character(0), rowPosition = character(0),
                     colPosition = character(0), residual = numeric(0), similar = integer(0))
  nr <- nrow(TD)
  indicator <- TD[, traits, drop = FALSE]
  similar <- rep(x = 2, times = nr)
  if (verbose) {
    cat("the residual method are large standardized residuals reported by the mixed
        model analysis\n")
  }
  for (ii in 1:length(traits)) {
    if (is.null(rLimit)) {
      if (rDf[ii] <= 20) {
        rLimit <- 2
      } else if (rDf[ii] <= 15773) {
        rLimit <- qnorm(p = 1 - 0.5 / rDf)
      } else {
        rLimit <- 4
      }
    }
    # identifies outliers as large standardized residuals
    tInd <- indicator[, traits[ii]] <- absResids[, traits[ii]] > rLimit
    # list all factors as similar to those with the large residuals
    nCFac <- length(commonFactor)
    genoRLarge <- ttInd <- rep(x = FALSE, times = nr)
    for (jj in 1:nCFac) {
      ttInd <- TD[tInd, commonFactor[jj]]
      ttInd <- droplevels(ttInd)
      ttInd <- levels(ttInd)
      ttInd <- TD[[commonFactor[jj]]] %in% ttInd
      genoRLarge <- genoRLarge + ttInd
    }
    genoRLarge <- genoRLarge > 0
    similar[genoRLarge] <- 1
    similar[tInd] <- 0
    nn <- sum(genoRLarge)
    if (nn > 0) {
      pMat0 <- data.frame(traitName = character(nn), traitValue = numeric(nn),
                          genotype = character(nn), entry = character(nn),
                          plotNo = character(nn), replicate = character(nn),
                          subBlock = character(nn),
                          rowId = character(nn), column = character(nn),
                          rowPosition = character(nn), colPosition = character(nn),
                          residual = numeric(nn), similar = integer(nn))
      pMat0$traitValue <- TD[genoRLarge, traits[ii]]
      pMat0$traitName <- rep(x = traits[ii], times = nn)
      pMat0$genotype <- TD[genoRLarge, "genotype"]
      pMat0$residual <- stdRes[genoRLarge, ii]
      pMat0$similar <- similar[genoRLarge]
      if (!is.null(entry)) {
        pMat0$entry <- TD[genoRLarge, entry]
      } else {
        pMat0$entry <- rep(NA, nn)
      }
      if (!is.null(plotNo)) {
        pMat0$plotNo <- TD[genoRLarge, plotNo]
      } else {
        pMat0$plotNo <- rep(NA, nn)
      }
      if (!is.null(repId)) {
        pMat0$replicate <- TD[genoRLarge, repId]
      } else {
        pMat0$replicate <- rep(NA, nn)
      }
      if (!is.null(subBlock)) {
        pMat0$subBlock <- TD[genoRLarge, subBlock]
      } else {
        pMat0$subBlock <- rep(NA, nn)
      }
      if (!is.null(rowId)) {
        pMat0$rowId <- TD[genoRLarge, rowId]
      } else {
        pMat0$rowId <- rep(NA, nn)
      }
      if (!is.null(colId)) {
        pMat0$column <- TD[genoRLarge, colId]
      } else {
        pMat0$column <- rep(NA, nn)
      }
      if (!is.null(rowCoordinates)) {
        pMat0$rowPosition <- TD[genoRLarge, rowCoordinates]
      } else {
        pMat0$rowPosition <- rep(NA, nn)
      }
      if (!is.null(colCoordinates)) {
        pMat0$colPosition <- TD[genoRLarge, colCoordinates]
      } else {
        pMat0$colPosition <- rep(NA, nn)
      }
      pMat <- rbind(pMat, pMat0)
    }
  }
  if (nrow(pMat) > 0 && verbose) {
  #  naNames0 <- c(FALSE, FALSE, FALSE, naNames, FALSE, FALSE)
  #  print(format(pMat[, !naNames0]), quote = FALSE)
  }
  return(indicator)
}
