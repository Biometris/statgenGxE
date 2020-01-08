## ----setup, include = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.dim = c(7, 4)
)
library(statgenGxE)
## Call requireNamespace here to prevent license output in first call in vignette.
requireNamespace("asreml", quietly = TRUE)
options(width = 90) 

## ----loadData---------------------------------------------------------------------------
data(dropsPheno)

## ----createTD---------------------------------------------------------------------------
## Create a TD object from dropsPheno
dropsTD <- createTD(data = dropsPheno, genotype = "Variety_ID", trial = "Experiment")

## ----TDbox------------------------------------------------------------------------------
## Create a box plot of dropsTD
## Color the boxes based on the variable scenarioFull.
## Plot in  descending order.
plot(dropsTD, plotType = "box", traits = "grain.yield", colorBy = "scenarioFull", 
     orderBy = "descending")

## ----geFW-------------------------------------------------------------------------------
## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "grain.yield")
summary(dropsFW)

## ----plotFWScatter, fig.width=5, fig.height=5, fig.show="hold"--------------------------
## Create scatter plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "scatter")

## ----plotFWLine, fig.width=5, fig.height=5, fig.show="hold"-----------------------------
## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line")

## ----plotFWTrellis, fig.width=5, fig.height=5, fig.show="hold"--------------------------
## Create trellis plot for Finlay Wilkinson analysis.
## Restrict to first 5 genotypes.
plot(dropsFW, plotType = "trellis", genotypes = c("11430", "A3", "A310", "A347", "A374"))

## ----geFWRep, eval=FALSE----------------------------------------------------------------
#  report(dropsFW, outfile = "./myReports/FWReport.pdf")

## ----geAmmi-----------------------------------------------------------------------------
## Run gxeAmmi with default settings.
dropsAm <- gxeAmmi(TD = dropsTD, trait = "grain.yield")
summary(dropsAm)

## ----geAmmi2----------------------------------------------------------------------------
## Run gxeAmmi. Let algorithm determine number of principal components.
dropsAm2 <- gxeAmmi(TD = dropsTD, trait = "grain.yield", nPC = NULL)
summary(dropsAm2)

## ----geAmmi3----------------------------------------------------------------------------
## Run gxeAmmi with three principal components.
## Exclude genotypes 11430 and A3.
dropsAm3 <- gxeAmmi(TD = dropsTD, trait = "grain.yield", nPC = 3, 
                    excludeGeno = c("11430", "A3"))

## ----geAmmiYear-------------------------------------------------------------------------
## Run gxeAmmi per year in the data.
dropsAmYear <- gxeAmmi(TD = dropsTD, trait = "grain.yield", byYear = TRUE)

## ----plotAmmi1, fig.width=5, fig.height=5, out.width="75%"------------------------------
## Create an AMMI1 and AMMI2 biplot.
plot(dropsAm, scale = 0.5, plotType = "AMMI1")

## ----plotAmmi2, fig.width=5, fig.height=5, out.width="75%"------------------------------
## Create an AMMI1 and AMMI2 biplot.
plot(dropsAm, scale = 0.5, plotType = "AMMI2")

## ----plotAmmiCol, fig.width=5, fig.height=5, out.width="75%"----------------------------
## Create an AMMI2 biplot.
## Color genotypes based on variable genetic_group. Use custom colors.
## Color environments base on variable scenarioFull
plot(dropsAm, scale = 0.4, plotType = "AMMI2", 
     colorGenoBy = "genetic_group", colGeno = c("red", "blue", "green", "yellow"),
     colorEnvBy = "scenarioFull")


## ----plotAmmiConvHull, fig.width=5, fig.height=5, out.width="75%"-----------------------
## Create an AMMI2 biplot with convex hull around the genotypes.
plot(dropsAm, scale = 0.4, plotType = "AMMI2", plotConvHull = TRUE, colorEnvBy = "scenarioFull")

## ----geAMMIRep, eval=FALSE--------------------------------------------------------------
#  report(dropsAm, outfile = "./myReports/AMMIReport.pdf")

## ----geGGE------------------------------------------------------------------------------
## Run gxeAmmi with default settings.
dropsGGE <- gxeGGE(TD = dropsTD, trait = "grain.yield")
summary(dropsGGE)

## ----plotGGE, fig.width=5, fig.height=5, out.width="75%"--------------------------------
## Create a GGE2 biplot.
plot(dropsGGE, scale = 0.5, plotType = "GGE2", plotConvHull = TRUE)

## ----geMegaEnv, R.options=list(digits=3)------------------------------------------------
## Compute mega environments.
dropsMegaEnv <- gxeMegaEnv(TD = dropsTD, trait = "grain.yield")

## ----geMegaEnvPred, R.options=list(digits=3)--------------------------------------------
if (requireNamespace(package = "asreml", quietly = TRUE)) {
  ## Compute BLUPs.
  ## Use asreml as engine for fitting model.
  geMegaEnvPred <- gxeTable(TD = dropsMegaEnv, trait = "grain.yield", engine = "asreml")
  ## Display BLUPs and associated standard errors.
  print(head(geMegaEnvPred$predictedValue))
  print(head(geMegaEnvPred$standardError))
}

## ----geStab, R.options=list(digits=3)---------------------------------------------------
## Compute stability measures for dropsTD.
dropsStab <- gxeStability(TD = dropsTD, trait = "grain.yield")
## In the summary print the top two percent of the genotypes.
summary(dropsStab, pctGeno = 2)

## ----plotStab---------------------------------------------------------------------------
plot(dropsStab)

## ----geStabRep, eval=FALSE--------------------------------------------------------------
#  report(dropsStab, outfile = "./myReports/stabReport.pdf")

## ----geStabMegaEnv, R.options=list(digits=3)--------------------------------------------
## Compute stability measures based on mega environments computed in the 
## previous paragraph.
dropsStabME <- gxeStability(TD = dropsMegaEnv, trait = "grain.yield", useMegaEnv = TRUE)
summary(dropsStabME, pctGeno = 2)

## ----geVClme----------------------------------------------------------------------------
## Use lme4 for fitting the models - only compound symmetry.
dropsVC <- gxeVarComp(TD = dropsTD, trait = "grain.yield")
summary(dropsVC)

## ----geVCasreml-------------------------------------------------------------------------
## Use asreml for fitting the models - eight models fitted. 
## Use AIC as criterion for determining the best model.
if (requireNamespace("asreml", quietly = TRUE)) {
  dropsVC2 <- gxeVarComp(TD = dropsTD, trait = "grain.yield", engine = "asreml",
                         criterion = "AIC")
  summary(dropsVC2)
}

## ----geVCPlot---------------------------------------------------------------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  plot(dropsVC2)
}

## ----geVCRep, eval=FALSE----------------------------------------------------------------
#  report(dropsVC2, outfile = "./myReports/varCompReport.pdf")

