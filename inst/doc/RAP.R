## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>"
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
wheatMeta$trPlotLength = 2
wheatMeta$trPlotWidth = 1
wheatTD <- setMeta(TD = wheatTD, meta = wheatMeta)

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

