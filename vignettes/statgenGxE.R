## ----setup, include = FALSE-------------------------------------------------------------
knitr::opts_chunk$set(
collapse = TRUE,
comment = "#>",
fig.dim = c(7, 4)
)
library(statgenGxE)
op <- options(width = 90)
options("statgen.trialColors" = c("#9E0142FF", "#35B779FF", "#B4DE2CFF",
                                  "#006837FF", "#D53E4FFF"))

## ----loadData---------------------------------------------------------------------------
data(dropsPheno)

## ----createTD---------------------------------------------------------------------------
## Create a TD object from dropsPheno.
dropsTD <- statgenSTA::createTD(data = dropsPheno, genotype = "Variety_ID", trial = "Experiment")

## ----TDbox------------------------------------------------------------------------------
## Create a box plot of dropsTD.
## Color the boxes based on the variable scenarioFull.
## Plot in  descending order.
plot(dropsTD, plotType = "box", traits = "grain.yield", colorTrialBy = "scenarioFull",
     orderBy = "descending")

## ----TDscatter, fig.dim = c(8.5, 8.5)---------------------------------------------------
## Create a scatter plot of dropsTD.
## Color the genotypes based on the variable geneticGroup.
## Color the histograms for trials based on the variable scenarioFull.
plot(dropsTD, plotType = "scatter", traits = "grain.yield", colorGenoBy = "geneticGroup", 
     colorTrialBy = "scenarioFull", 
     trialOrder = c("Gai12W", "Kar13R", "Kar12W", "Kar13W", "Mar13R", "Mur13W",
                    "Mur13R", "Ner12R", "Cam12R", "Cra12R"))

## ----colorOpts, eval=FALSE--------------------------------------------------------------
#  ## Set default colors for genotypes and trials.
#  options("statgen.genoColors" = c("blue", "green", "yellow"))
#  options("statgen.trialColors" = c("red", "brown", "purple"))

## ----geVarComp, message=FALSE-----------------------------------------------------------
## Fit a model where trials are nested within scenarios.
dropsVarComp <- gxeVarComp(TD = dropsTD, trait = "grain.yield", nestingFactor = "scenarioFull")
summary(dropsVarComp)

## ----diag, eval=FALSE-------------------------------------------------------------------
#  ## Print diagnostics - output suppressed because of the large number of rows.
#  diagnostics(dropsVarComp)

## ----vcHerit, R.options=list(digits=4)--------------------------------------------------
## Extract variance components.
vc(dropsVarComp)
## Compute heritability.
herit(dropsVarComp)

## ----VarCompPlot------------------------------------------------------------------------
## Plot the results of the fitted model.
plot(dropsVarComp)

## ----predict, R.options=list(digits=4)--------------------------------------------------
## Predictions of the genotype main effect.
predGeno <- predict(dropsVarComp)
head(predGeno)
## predictions at the level of genotype x scenarioFull.
predGenoTrial <- predict(dropsVarComp, predictLevel = "scenarioFull")
head(predGenoTrial)

## ----geFW, R.options=list(digits=6)-----------------------------------------------------
## Perform a Finlay-Wilkinson analysis for all trials.
dropsFW <- gxeFw(TD = dropsTD, trait = "grain.yield")
summary(dropsFW)

## ----plotFWScatter, fig.width=5, fig.height=4-------------------------------------------
## Create scatter plot for Finlay Wilkinson analysis.
## Color genotypes by geneticGroup.
plot(dropsFW, plotType = "scatter", colorGenoBy = "geneticGroup")

## ----plotFWLine, fig.width=5, fig.height=4----------------------------------------------
## Create line plot for Finlay Wilkinson analysis.
## Color genotypes by geneticGroup.
plot(dropsFW, plotType = "line", colorGenoBy = "geneticGroup")

## ----plotFWTrellis, fig.width=5, fig.height=4-------------------------------------------
## Create trellis plot for Finlay Wilkinson analysis.
## Restrict to first 5 genotypes.
plot(dropsFW, plotType = "trellis", genotypes = c("11430", "A3", "A310", "A347", "A374"))

## ----plotFWScatterFit, fig.width=5, fig.height=4----------------------------------------
## Create scatter plot of fitted values for Finlay Wilkinson analysis.
## Color genotypes by geneticGroup.
plot(dropsFW, plotType = "scatterFit", colorGenoBy = "geneticGroup")

## ----geAmmi, R.options=list(digits=6)---------------------------------------------------
## Run gxeAmmi for grain.yield.
dropsAm <- gxeAmmi(TD = dropsTD, trait = "grain.yield")
summary(dropsAm)

## ----geAmmi2, R.options=list(digits=6)--------------------------------------------------
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
## Create an AMMI1 plot.
plot(dropsAm, plotType = "AMMI1")

## ----plotAmmi2, fig.width=5, fig.height=5, out.width="75%"------------------------------
## Create an AMMI2 biplot with symmetric scaling.
plot(dropsAm, scale = 0.5, plotType = "AMMI2")

## ----plotAmmiCol, fig.width=5, fig.height=5, out.width="75%"----------------------------
## Create an AMMI2 biplot.
## Color genotypes based on variable geneticGroup. Use custom colors.
## Color environments based on variable scenarioFull
plot(dropsAm, scale = 0.4, plotType = "AMMI2", 
     colorGenoBy = "geneticGroup", colGeno = c("red", "blue", "green", "yellow"),
     colorEnvBy = "scenarioFull")


## ----plotAmmiConvHull, fig.width=5, fig.height=5, out.width="75%"-----------------------
## Create an AMMI2 biplot with convex hull around the genotypes.
plot(dropsAm, scale = 0.4, plotType = "AMMI2", plotConvHull = TRUE, colorEnvBy = "scenarioFull")

## ----geGGE, R.options=list(digits=6)----------------------------------------------------
## Run gxeGGE with default settings.
dropsGGE <- gxeGGE(TD = dropsTD, trait = "grain.yield")
summary(dropsGGE) 

## ----plotGGE, fig.width=5, fig.height=5, out.width="75%"--------------------------------
## Create a GGE2 biplot.
plot(dropsGGE, plotType = "GGE2")

## ----geMegaEnv, R.options=list(digits=5)------------------------------------------------
## Compute mega environments.
dropsMegaEnv <- gxeMegaEnv(TD = dropsTD, trait = "grain.yield")
## Summarize results.
summary(dropsMegaEnv)

## ----geMegaEnvPred, R.options=list(digits=5)--------------------------------------------
if (requireNamespace(package = "asreml", quietly = TRUE)) {
  ## Compute BLUPs.
  ## Use asreml as engine for fitting model.
  geMegaEnvPred <- predict(dropsMegaEnv, engine = "asreml")
  ## Display BLUPs.
  head(geMegaEnvPred$predictedValue)
}
if (requireNamespace(package = "asreml", quietly = TRUE)) {
  ## Display standard errors of the BLUPs.
  head(geMegaEnvPred$standardError)
}

## ----scatterMegaEnv---------------------------------------------------------------------
if (requireNamespace(package = "asreml", quietly = TRUE)) {
  ## Create a scatter plot of predictions in mega environments.
  ## Color genotypes based on geneticGroup.
  plot(dropsMegaEnv, engine = "asreml", colorGenoBy = "geneticGroup")
}

## ----geStab, R.options=list(digits=6)---------------------------------------------------
## Compute stability measures for dropsTD.
dropsStab <- gxeStability(TD = dropsTD, trait = "grain.yield")
## In the summary print the top two percent of the genotypes.
summary(dropsStab, pctGeno = 2)

## ----plotStab---------------------------------------------------------------------------
## Create plots for the different stability measures.
## Color genotypes by geneticGroup.
plot(dropsStab, colorGenoBy = "geneticGroup")

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

## ----fitRes-----------------------------------------------------------------------------
## Extract the fitted values and residuals for the Finlay-Wilkinson model.
fitFW <- fitted(dropsFW)
resFW <- residuals(dropsFW)

## Create a diagnostic plot of fitted values against residuals.
fitResFW <- merge(fitFW, resFW, by = c("trial", "genotype"))
ggplot2::ggplot(fitResFW, ggplot2::aes(x = fittedValue, y = residual,
                                       color = trial)) +
  ggplot2::geom_point()


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

## ----winddown, include = FALSE------------------------------------------------
options(op)

