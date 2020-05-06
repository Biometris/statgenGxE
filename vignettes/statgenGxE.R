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
## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = dropsPheno, genotype = "Variety_ID", trial = "Experiment")

## ----TDbox------------------------------------------------------------------------------
## Create a box plot of dropsTD.
## Color the boxes based on the variable scenarioFull.
## Plot in  descending order.
plot(dropsTD, plotType = "box", traits = "grain.yield", colorBy = "scenarioFull", 
     orderBy = "descending")

## ----TDscatter, fig.dim = c(7, 7)-------------------------------------------------------
## Create a scatter plot of dropsTD.
## Color the genotypes based on the variable genetic_group.
plot(dropsTD, plotType = "scatter", traits = "grain.yield", colorBy = "genetic_group")

## ----geFW-------------------------------------------------------------------------------
## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "grain.yield")
summary(dropsFW)

## ----plotFWScatter, fig.width=5, fig.height=4-------------------------------------------
## Create scatter plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "scatter")

## ----plotFWLine, fig.width=5, fig.height=4----------------------------------------------
## Create line plot for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "line")

## ----plotFWTrellis, fig.width=5, fig.height=4-------------------------------------------
## Create trellis plot for Finlay Wilkinson analysis.
## Restrict to first 5 genotypes.
plot(dropsFW, plotType = "trellis", genotypes = c("11430", "A3", "A310", "A347", "A374"))

## ----plotFWScatterFit, fig.width=5, fig.height=4----------------------------------------
## Create scatter plot of fitted values for Finlay Wilkinson analysis.
plot(dropsFW, plotType = "scatterFit")

## ----geAmmi-----------------------------------------------------------------------------
## Run gxeAmmi for grain.yield.
## Scale the residuals before running the principal component analysis.
dropsAm <- gxeAmmi(TD = dropsTD, trait = "grain.yield", scale = TRUE)
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
## Create an AMMI1 biplot.
plot(dropsAm, plotType = "AMMI1")

## ----plotAmmi2, fig.width=5, fig.height=5, out.width="75%"------------------------------
## Create an AMMI2 biplot.
plot(dropsAm, scale = 0.5, plotType = "AMMI2")

## ----plotAmmiCol, fig.width=5, fig.height=5, out.width="75%"----------------------------
## Create an AMMI2 biplot.
## Color genotypes based on variable genetic_group. Use custom colors.
## Color environments based on variable scenarioFull
plot(dropsAm, scale = 0.4, plotType = "AMMI2", 
     colorGenoBy = "genetic_group", colGeno = c("red", "blue", "green", "yellow"),
     colorEnvBy = "scenarioFull")


## ----plotAmmiConvHull, fig.width=5, fig.height=5, out.width="75%"-----------------------
## Create an AMMI2 biplot with convex hull around the genotypes.
plot(dropsAm, scale = 0.4, plotType = "AMMI2", plotConvHull = TRUE, colorEnvBy = "scenarioFull")

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
## Summarize results.
summary(dropsMegaEnv)

## ----geMegaEnvPred, R.options=list(digits=3)--------------------------------------------
if (requireNamespace(package = "asreml", quietly = TRUE)) {
  ## Compute BLUPs.
  ## Use asreml as engine for fitting model.
  geMegaEnvPred <- predict(dropsMegaEnv, engine = "asreml")
  ## Display BLUPs and associated standard errors.
  print(head(geMegaEnvPred$predictedValue))
  print(head(geMegaEnvPred$standardError))
}

## ----scatterMegaEnv---------------------------------------------------------------------
if (requireNamespace(package = "asreml", quietly = TRUE)) {
  ## Create a scatter plot of predictions in mega environments.
  ## Color genotypes based on genetic_group.
  plot(dropsMegaEnv, engine = "asreml", colorBy = "genetic_group")
}

## ----geStab, R.options=list(digits=3)---------------------------------------------------
## Compute stability measures for dropsTD.
dropsStab <- gxeStability(TD = dropsTD, trait = "grain.yield")
## In the summary print the top two percent of the genotypes.
summary(dropsStab, pctGeno = 2)

## ----plotStab---------------------------------------------------------------------------
plot(dropsStab)

## ----geStabMegaEnv, R.options=list(digits=3), eval = FALSE------------------------------
#  ## Compute stability measures based on mega environments computed in the
#  ## previous paragraph.
#  dropsStabME <- gxeStability(TD = dropsMegaEnv, trait = "grain.yield", useMegaEnv = TRUE)
#  summary(dropsStabME, pctGeno = 2)

## ----geVClme----------------------------------------------------------------------------
## Use lme4 for fitting the models - only compound symmetry.
dropsVC <- gxeVarCov(TD = dropsTD, trait = "grain.yield")
summary(dropsVC)

## ----geVCasreml-------------------------------------------------------------------------
## Use asreml for fitting the models - eight models fitted. 
## Use AIC as criterion for determining the best model.
if (requireNamespace("asreml", quietly = TRUE)) {
  dropsVC2 <- gxeVarCov(TD = dropsTD, trait = "grain.yield", engine = "asreml", criterion = "AIC")
  summary(dropsVC2)
}

## ----geVCPlot---------------------------------------------------------------------------
if (requireNamespace("asreml", quietly = TRUE)) {
  plot(dropsVC2)
}

## ----reports, eval=FALSE----------------------------------------------------------------
#  ## Create a report for the Finlay Wilkinson analysis.
#  report(dropsFW, outfile = "./myReports/FWReport.pdf")
#  
#  ## Create a report for the AMMI analysis.
#  report(dropsAm, outfile = "./myReports/AMMIReport.pdf")
#  
#  ## Create a report for the GGE analysis.
#  report(dropsGGE, outfile = "./myReports/GGEReport.pdf")
#  
#  ## Create a report for the stability analysis.
#  report(dropsStab, outfile = "./myReports/stabReport.pdf")
#  
#  ## Create a report for the analysis of two-way GxE tables.
#  report(dropsVC2, outfile = "./myReports/varCompReport.pdf")

