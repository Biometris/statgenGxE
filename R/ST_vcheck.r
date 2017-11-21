#' Identifying outliers
#'
#' Identifies large standardized residuals from a REML analysis.
#'
#' @inheritParams ST.run.model
#'
#' @param stdRes A data frame object containing the standardised residuals obtained
#' from the mixed modelling analysis.
#' @param rDf An integer (vector) specifying the residual degrees of freedom.
#' @param rLimit An integer specifying a limit for detection of large standardized residuals;
#' if this is not set, the limit is set automatically according to the number of residual
#' degrees of freedom.
#' @param entry A string specifying the column name of the entry numbers.
#' @param plotNo A string specifying the column name of the plot numbers.
#' @param commonFactor A string vector specifying factors to define similar units; default,
#' \code{commonFactor = genotype}.
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
#'                          rep = "Rep", row = "Row", col = "Column",
#'                          engine = "lme4") #engine = "asreml"
#' stdResid <- resid(myModel$mFix, scaled = TRUE)
#' rDf <- nrow(model.frame(myModel$mFix)) - length(lme4::fixef(myModel$mFix))
#' vCheck <- ST.vcheck(TD = myTD, stdRes = stdResid, rDf = rDf, trait = "yield",
#'                     rep = "Rep", plotNo = "Plot", row = "Row", col = "Column",
#'                     verbose = TRUE, commonFactor = c("genotype", "Row"))
#'
#' @export
ST.vcheck <- function(TD,
                      stdRes,
                      rDf,
                      rLimit = NA,
                      trait,
                      entry = NA,
                      plotNo = NA,
                      rep = NA,
                      subBlock = NA,
                      row = NA,
                      col = NA,
                      rowCoordinates = NA,
                      colCoordinates = NA,
                      commonFactor = "genotype",
                      verbose = FALSE) {
  nr <- nrow(TD)
  if (is.vector(stdRes)) {
    stdRes <- data.frame(stdRes)
    names(stdRes) <- trait
  }
  absResids <- abs(stdRes)
  pMat <- data.frame(traitName = character(0), traitValue = numeric(0),
                     genotype = character(0), entry = character(0), plotNo = character(0),
                     replicate = character(0), subBlock = character(0), row = character(0),
                     column = character(0), rowPosition = character(0),
                     colPosition = character(0), residual = numeric(0), similar = integer(0))
  testNames <- c(entry, plotNo, rep, subBlock, row, col, rowCoordinates, colCoordinates)
  naNames <- is.na(testNames)
  testNames <- testNames[!naNames]
  if (!all(testNames %in% names(TD))) {
    stop(paste(testNames[!(testNames %in% names(TD))], collapse = ","),
         " not found in the names of data.\n")
  }
  indicator <- TD[, trait, drop = FALSE]
  similar <- rep(x = 2, times = nr)
  if (verbose) {
    cat("the residual method are large standardized residuals reported by the mixed
        model analysis\n")
  }
  if (length(trait) != length(rDf)) {
    stop("number of traits does not equal to number of residual degrees of freedom.\n")
  }
  if (ncol(absResids) != length(rDf)) {
    stop("number of columns in standardised residuals does not equal to number of residual
         degrees of freedom.\n")
  }
  if (ncol(absResids) != length(trait)) {
    stop("number of traits does not equal to number of columns in standardised residuals.\n")
  }
  for (ii in 1:length(trait)) {
    if (is.na(rLimit)) {
      if (rDf[ii] <= 20) {
        rLimit <- 2
      } else if (rDf[ii] <= 15773) {
        rLimit <- qnorm(p = 1 - 0.5 / rDf)
      } else {
        rLimit <- 4
      }
    }
    # identifies outliers as large standardized residuals
    tInd <- indicator[, trait[ii]] <- absResids[, trait[ii]] > rLimit
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
                          row = character(nn), column = character(nn),
                          rowPosition = character(nn), colPosition = character(nn),
                          residual = numeric(nn), similar = integer(nn))
      pMat0$traitValue <- TD[genoRLarge, trait[ii]]
      pMat0$traitName <- rep(x = trait[ii], times = nn)
      pMat0$genotype <- TD[genoRLarge, "genotype"]
      pMat0$residual <- stdRes[genoRLarge, ii]
      pMat0$similar <- similar[genoRLarge]
      if (!is.na(entry)) {
        pMat0$entry <- TD[genoRLarge, entry]
      } else {
        pMat0$entry <- rep(NA, nn)
      }
      if (!is.na(plotNo)) {
        pMat0$plotNo <- TD[genoRLarge, plotNo]
      } else {
        pMat0$plotNo <- rep(NA, nn)
      }
      if (!is.na(rep)) {
        pMat0$replicate <- TD[genoRLarge, rep]
      } else {
        pMat0$replicate <- rep(NA, nn)
      }
      if (!is.na(subBlock)) {
        pMat0$subBlock <- TD[genoRLarge, subBlock]
      } else {
        pMat0$subBlock <- rep(NA, nn)
      }
      if (!is.na(row)) {
        pMat0$row <- TD[genoRLarge, row]
      } else {
        pMat0$row <- rep(NA, nn)
      }
      if (!is.na(col)) {
        pMat0$column <- TD[genoRLarge, col]
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
    naNames0 <- c(FALSE, FALSE, FALSE, naNames, FALSE, FALSE)
    print(format(pMat[, !naNames0]), quote = FALSE)
  }
  return(indicator)
}
