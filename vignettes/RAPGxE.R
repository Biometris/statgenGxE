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
wheatTD <- createTD(data = wheatChl, 
                    genotype = "trt", repId = "rep", subBlock = "bl", 
                    rowCoord = "row", colCoord = "col")
modWheatSpTot <- STRunModel(TD = wheatTD, traits = "GY", what = "fixed", design = "res.rowcol")
## Create a TD object containing BLUEs and standard errors of BLUEs.
TDGxE <- SSAtoTD(SSA = modWheatSpTot, what = c("BLUEs", "seBLUEs"))

## ----geVClme-------------------------------------------------------------
## Use lme4 for fitting the models - only compound symmetry
geVC <- gxeVarComp(TD = TDGxE, trait = "BLUEs_GY")
summary(geVC)

## ----geVCasreml----------------------------------------------------------
## Use asreml for fitting the models - 8 models fitted. 
## Use AIC as criterion for determining the best model.
if (requireNamespace("asreml")) {
  geVC2 <- gxeVarComp(TD = TDGxE, trait = "BLUEs_GY", engine = "asreml", 
                      criterion = "AIC")
  summary(geVC2)
}

## ----geVCPlot, out.width="75%"-------------------------------------------
if (requireNamespace("asreml")) {
  plot(geVC2)
}

## ----geVCRep, eval=FALSE-------------------------------------------------
#  report(geVC2, outfile = "./myReports/varCompReport.pdf")

## ----geAmmi--------------------------------------------------------------
## Run gxeAmmi with default settings.
geAm <- gxeAmmi(TD = TDGxE, trait = "BLUEs_GY")
summary(geAm)

## ----geAmmi2-------------------------------------------------------------
## Run gxeAmmi with 3 principal components.
geAm2 <- gxeAmmi(TD = TDGxE, trait = "BLUEs_GY", nPC = 3)
summary(geAm2)

## ----geAmmi3-------------------------------------------------------------
## Run gxeAmmi. Algorithm determines number of principal components.
geAm3 <- gxeAmmi(TD = TDGxE, trait = "BLUEs_GY", nPC = NULL)
summary(geAm3)

## ----plotAmmi,fig.width=5,fig.height=5,out.width="47%",fig.show="hold"----
## Create an AMMI1 and AMMI2 biplot.
plot(geAm, scale = 0.5, plotType = "AMMI1")
plot(geAm, scale = 0.5, plotType = "AMMI2")

## ----geAMMIRep, eval=FALSE-----------------------------------------------
#  report(geAm, outfile = "./myReports/AMMIReport.pdf")

## ----geFW----------------------------------------------------------------
geFW <- gxeFw(TD = TDGxE, trait = "BLUEs_GY")
summary(geFW)

## ----plotFW,fig.width=5,fig.height=5,out.width="47%",fig.show="hold"-----
plot(geFW, plotType = "scatter")
plot(geFW, plotType = "line")
plot(geFW, plotType = "trellis")

## ----geFWRep, eval=FALSE-------------------------------------------------
#  report(geFW, outfile = "./myReports/FWReport.pdf")

## ----geMegaEnv-----------------------------------------------------------
geMegaEnv <- gxeMegaEnv(TD = TDGxE, trait = "BLUEs_GY")

## ----geStab--------------------------------------------------------------
geStab <- gxeStability(TD = TDGxE, trait = "BLUEs_GY")
summary(geStab)

## ----plotStab------------------------------------------------------------
plot(geStab)

## ----geStabRep, eval=FALSE-----------------------------------------------
#  report(geStab, outfile = "./myReports/stabReport.pdf")

## ----geStabMegaEnv-------------------------------------------------------
## Compute stabilities measures based on mega environments computed in 
## previous paragraph.
geStabME <- gxeStability(TD = geMegaEnv, trait = "BLUEs_GY", useMegaEnv = TRUE)
summary(geStabME)

