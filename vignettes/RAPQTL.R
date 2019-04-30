## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.dim = c(7, 5)
)
library(RAP)

## ---- include = FALSE----------------------------------------------------
## Recreate data from last step in RAP vignette.
data("wheatChl")
wheatTD <- createTD(data = wheatChl, genotype = "trt", repId = "rep", 
                    subBlock = "bl", rowCoord = "row", colCoord = "col")
modWheatSpTot <- fitTD(TD = wheatTD, traits = "GY", what = "fixed", 
                       design = "res.rowcol")

## ----createCross---------------------------------------------------------
cr1 <- SSAtoCross(SSA = modWheatSpTot, trial = "C_SWS_12",  
                  what = "BLUEs", genoFile = system.file("extdata", "markersWheatChl.csv", package = "RAP"), na.strings = c("N", "H"))

## ----qualContr-----------------------------------------------------------
#QTLMapQC(cross = cr1)

