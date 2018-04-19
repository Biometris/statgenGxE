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
wheatMeta <- getMeta(TD = wheatTD)
head(wheatMeta)

## ----setMeta-------------------------------------------------------------
wheatMeta$trLocation <- "Santa Rosa"
wheatMeta$trDate <- as.Date(rep(c("310811", "310812"), times = 2), "%d%m%y")
wheatMeta$trLat <- -36.32
wheatMeta$trLong <- -71.55
wheatMeta$trPlotWidth = 2
wheatMeta$trPlotLength = 1
wheatTD <- setMeta(TD = wheatTD, meta = wheatMeta)

## ----addTD---------------------------------------------------------------
wheatTD <- addTD(TD = wheatTD, data = wheatChl[wheatChl$trial == "C_SWS_12", ], 
                 genotype = "trt_id", repId = "rep", subBlock = "bl", 
                 rowId = "row", colId = "col", rowCoord = "row", 
                 colCoord = "col", trLocation = "Cauquenes", 
                 trDate = as.Date("070712", "%d%m%y"), trLat = -35.58,
                 trLong = -72.17, trPlotWidth = 2, trPlotLength = 1)
## Inspect the meta data after the extra trial was added.
getMeta(TD = wheatTD)

## ----TDsum---------------------------------------------------------------
summary(wheatTD, trial = "SR_FI_11", traits = "GY")

## ----layoutPlot----------------------------------------------------------
plot(wheatTD, trials = "SR_FI_11")

## ----mapPlot-------------------------------------------------------------
plot(wheatTD, plotType = "map")

## ----fitSp, message=FALSE------------------------------------------------
modWheat <- STRunModel(TD = wheatTD, trials = "SR_FI_11", traits = "GY",
                       design = "res.rowcol")

## ----fitSum, message=FALSE-----------------------------------------------
## Set nBest to 5 to decrease size of output.
summary(modWheat, nBest = 5)

## ----basePlot------------------------------------------------------------
plot(modWheat, what = "random")

## ----spatPlot------------------------------------------------------------
plot(modWheat, plotType = "spatial")

## ----modRep, eval=FALSE--------------------------------------------------
#  ## Create a report in the current working directory
#  report(modWheat)
#  ## Create a report for the model with genotype fitted as random.
#  report(modWheat, outfile = "./myReports/wheatReport.pdf", what = "random")

## ----extBLUEs------------------------------------------------------------
## Extract BLUEs
BLUEsWheat <- STExtract(SSA = modWheat, what = "BLUEs")
## Extract BLUEs and BLUPs
predWheat <- STExtract(SSA = modWheat, what = c("BLUEs", "BLUPs"))

## ----extBLUEsKeep--------------------------------------------------------
BLUEsWheat2 <- STExtract(SSA = modWheat, what = "BLUEs", keep = "trial")

