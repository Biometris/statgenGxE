## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.dim = c(7, 5)
)
library(RAP)

## ----createTD------------------------------------------------------------
data("wheatChl")
wheatTD <- createTD(data = wheatChl[wheatChl$trial != "C_SWS_12", ], 
                    genotype = "trt_id", repId = "rep", subBlock = "bl", 
                    rowId = "row", colId = "col", rowCoord = "row", 
                    colCoord = "col")

## ----getMeta-------------------------------------------------------------
(wheatMeta <- getMeta(TD = wheatTD))

## ----setMeta-------------------------------------------------------------
wheatMeta$trLocation <- "Santa Rosa"
wheatMeta$trDate <- as.Date(rep(c("310811", "310812"), times = 2), "%d%m%y")
wheatMeta$trLat <- -36.32
wheatMeta$trLong <- -71.55
wheatMeta$trPlWidth = 2
wheatMeta$trPlLength = 1
wheatTD <- setMeta(TD = wheatTD, meta = wheatMeta)

## ----addTD, R.options=list(width=90)----------------------------------------------------
wheatTD <- addTD(TD = wheatTD, data = wheatChl[wheatChl$trial == "C_SWS_12", ], 
                 genotype = "trt_id", repId = "rep", subBlock = "bl", 
                 rowId = "row", colId = "col", rowCoord = "row", 
                 colCoord = "col", trLocation = "Cauquenes", 
                 trDate = as.Date("070912", "%d%m%y"), trLat = -35.58,
                 trLong = -72.17, trPlWidth = 2, trPlLength = 1)
## Inspect the meta data after the extra trial was added.
getMeta(TD = wheatTD)

## ----TDsum---------------------------------------------------------------
summary(wheatTD, trial = "SR_FI_11", traits = "GY")

## ----layoutPlot----------------------------------------------------------
plot(wheatTD, trials = "SR_FI_11")

## ----mapPlot-------------------------------------------------------------
plot(wheatTD, plotType = "map")

## ----fitSp, message=FALSE------------------------------------------------
modWheatSp <- STRunModel(TD = wheatTD, trials = "SR_FI_11", traits = "GY",
                         design = "res.rowcol")

## ----fitSpSm, message=FALSE----------------------------------------------
## Fit a single trial model with genotype as random effect.
modWheatSp2 <- STRunModel(TD = wheatTD, trials = "SR_FI_11", traits = "GY",
                          what = "random", design = "res.rowcol")

## ----fitSpCtr, message=FALSE---------------------------------------------
modWheatSp3 <- STRunModel(TD = wheatTD, trials = "SR_FI_11", traits = "GY",
                          design = "res.rowcol", 
                          control = list(nSeg = c(20, 20)))

## ----fitAs, message=FALSE, results='hide'--------------------------------
if (requireNamespace("asreml")) {
  modWheatAs <- STRunModel(TD = wheatTD, trials = "SR_FI_11", traits = "GY",
                           design = "res.rowcol", trySpatial = TRUE,
                           engine = "asreml", control = list(criterion = "BIC"))
}

## ----spatCh--------------------------------------------------------------
if (requireNamespace("asreml")) {
  modWheatAs$SR_FI_11$spatial
}

## ----fitSum, message=FALSE-----------------------------------------------
## Set nBest to 5 to decrease size of output.
summary(modWheatSp, nBest = 5)

## ----basePlot------------------------------------------------------------
## Plot for the model with genotype fitted as random effect.
plot(modWheatSp, what = "random")

## ----spatPlot------------------------------------------------------------
## Spatial plot for the model with genotype fitted as fixed effect.
plot(modWheatSp, plotType = "spatial")

## ----modRep, eval=FALSE--------------------------------------------------
#  ## Create a report in the current working directory
#  report(modWheatSp)
#  ## Create a report for the model with genotype fitted as random.
#  report(modWheatSp, outfile = "./myReports/wheatReport.pdf", what = "random")

## ----extBLUEs------------------------------------------------------------
## Extract BLUEs
BLUEsWheat <- STExtract(SSA = modWheatSp, what = "BLUEs")
## Extract BLUEs and BLUPs
predWheat <- STExtract(SSA = modWheatSp, what = c("BLUEs", "BLUPs"))

## ----extBLUEsKeep--------------------------------------------------------
BLUEsWheat2 <- STExtract(SSA = modWheatSp, what = "BLUEs", keep = "trial")
head(BLUEsWheat2[[1]])

## ----extFit--------------------------------------------------------------
fitVals <- STExtract(SSA = modWheatSp, what = "fitted", 
                     keep = c("trial", "repId"))
head(fitVals[[1]])

