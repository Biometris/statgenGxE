## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.dim = c(6, 4)
)
library(statgenGxE)
## Call requireNamespace here to prevent license output in first call in vignette.
requireNamespace("asreml", quietly = TRUE)

## ---- include = FALSE, message = FALSE-----------------------------------
## Recreate data from last step in statgenSSA vignette.
data("wheatChl")
wheatTD <- createTD(data = wheatChl, genotype = "trt", repId = "rep", 
                    subBlock = "bl", rowCoord = "row", colCoord = "col")
modWheatSpTot <- statgenSSA::fitTD(TD = wheatTD, traits = "GY", what = "fixed", 
                                   design = "res.rowcol")
## Create a TD object containing BLUEs and standard errors of BLUEs.
TDGxE <- statgenSSA::SSAtoTD(SSA = modWheatSpTot, what = c("BLUEs", "seBLUEs"))

## ----geVClme-------------------------------------------------------------
## Use lme4 for fitting the models - only compound symmetry.
geVC <- gxeVarComp(TD = TDGxE, trait = "BLUEs_GY")
summary(geVC)

## ----geVCasreml----------------------------------------------------------
## Use asreml for fitting the models - eight models fitted. 
## Use AIC as criterion for determining the best model.
if (requireNamespace("asreml", quietly = TRUE)) {
  geVC2 <- gxeVarComp(TD = TDGxE, trait = "BLUEs_GY", engine = "asreml", 
                      criterion = "AIC")
  summary(geVC2)
}

## ----geVCPlot, out.width="75%"-------------------------------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  plot(geVC2)
}

## ----geVCRep, eval=FALSE-------------------------------------------------
#  report(geVC2, outfile = "./myReports/varCompReport.pdf")

## ----geAmmi--------------------------------------------------------------
## Run gxeAmmi with default settings.
geAm <- gxeAmmi(TD = TDGxE, trait = "BLUEs_GY")
summary(geAm)

## ----geAmmi2-------------------------------------------------------------
## Run gxeAmmi. Algorithm determines number of principal components.
geAm2 <- gxeAmmi(TD = TDGxE, trait = "BLUEs_GY", nPC = NULL)
summary(geAm2)

## ----geAmmi3-------------------------------------------------------------
## Run gxeAmmi with three principal components.
## Exclude genotypes G278 and G279.
geAm3 <- gxeAmmi(TD = TDGxE, trait = "BLUEs_GY", nPC = 3, 
                 excludeGeno = c("G278", "G279"))
summary(geAm3)

## ----plotAmmi, fig.width=5, fig.height=5, out.width="47%", fig.show="hold"----
## Create an AMMI1 and AMMI2 biplot.
plot(geAm, scale = 0.5, plotType = "AMMI1")
plot(geAm, scale = 0.5, plotType = "AMMI2")

## ----plotAmmi2, fig.width=5, fig.height=5, out.width="75%"---------------
## Create an AMMI2 biplot with convex hull around the genotypes and genotype names 
## displayed. Blow up genotypic scores by using envFactor = 0.3
plot(geAm, scale = 0.5, plotType = "AMMI2", sizeGeno = 2, plotConvHull = TRUE, 
     envFactor = 0.3)


## ----geAMMIRep, eval=FALSE-----------------------------------------------
#  report(geAm, outfile = "./myReports/AMMIReport.pdf")

## ----geGGE---------------------------------------------------------------
## Run gxeAmmi with default settings.
geGGE <- gxeAmmi(TD = TDGxE, trait = "BLUEs_GY", GGE = TRUE)
summary(geGGE)

## ----plotGGE, fig.width=5, fig.height=5, out.width="75%"-----------------
## Create an GGE1 and GGE2 biplot.
plot(geGGE, scale = 0.5, plotType = "GGE2", plotConvHull = TRUE)

## ----geFW----------------------------------------------------------------
## Perform a Finlay-Wilkinson analysis for all trials.
geFW <- gxeFw(TD = TDGxE, trait = "BLUEs_GY")
summary(geFW)

## ----plotFW,fig.width=5,fig.height=5,fig.show="hold"---------------------
plot(geFW, plotType = "scatter")
plot(geFW, plotType = "line")
plot(geFW, plotType = "trellis")

## ----geFWRep, eval=FALSE-------------------------------------------------
#  report(geFW, outfile = "./myReports/FWReport.pdf")

## ----geMegaEnv-----------------------------------------------------------
geMegaEnv <- gxeMegaEnv(TD = TDGxE, trait = "BLUEs_GY")

## ----geMegaEnvPred-------------------------------------------------------
if (requireNamespace(package = "asreml", quietly = TRUE)) {
  geMegaEnvPred <- gxeTable(TD = geMegaEnv, trait = "BLUEs_GY", engine = "asreml")
  head(geMegaEnvPred$predictedValue)
}

## ----geStab--------------------------------------------------------------
geStab <- gxeStability(TD = TDGxE, trait = "BLUEs_GY")
summary(geStab, pctGeno = 2)

## ----plotStab------------------------------------------------------------
plot(geStab)

## ----geStabRep, eval=FALSE-----------------------------------------------
#  report(geStab, outfile = "./myReports/stabReport.pdf")

## ----geStabMegaEnv-------------------------------------------------------
## Compute stabilities measures based on mega environments computed in the 
## previous paragraph.
geStabME <- gxeStability(TD = geMegaEnv, trait = "BLUEs_GY", useMegaEnv = TRUE)
summary(geStabME, pctGeno = 2)

