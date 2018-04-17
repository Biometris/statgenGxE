## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
library(RAP)

## ----createTD------------------------------------------------------------
data("wheatAus")
wheatTD <- createTD(data = wheatAus, genotype = "geno", trial = "Trial",
                    repId = "rep", rowId = "row", colId = "col", 
                    rowCoord = "row", colCoord = "col")

## ----getMeta-------------------------------------------------------------
wheatMeta <- getMeta(TD = wheatTD)
head(wheatMeta)

## ----setMeta-------------------------------------------------------------
wheatMeta$trDate <- rep(as.Date(c("01/03/08", "01/03/09"), "%d/%m/%y"), each = 7)
wheatMeta$trLat <- 33.8688
wheatMeta$trLong <- 151.2093
wheatTD <- setMeta(TD = wheatTD, meta = wheatMeta)

## ----TDsum---------------------------------------------------------------
summary(wheatTD, trial = "08CM", traits = "yield")
plot(wheatTD, trials = "08CM")

## ----mapPlot-------------------------------------------------------------
plot(wheatTD, plotType = "map")

