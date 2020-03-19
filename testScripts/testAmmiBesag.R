library(agridat)
library(statgenSTA)

## Compute BLUEs and wt for besag.met data from agridat.
## 64 genotypes - 6 environments.
besagTD <- createTD(besag.met, genotype = "gen", trial = "county", rowCoord = "row",
                    colCoord = "col", repId = "rep", subBlock = "block")
besagSTA <- fitTD(besagTD, traits = "yield", design = "rcbd")
besagBLUEs <- STAtoTD(besagSTA, what = c("BLUEs", "seBLUEs"), addWt = TRUE)

## Unweighted ammi with original function (imputation of missings + svd)
amOrig <- gxeAmmiOrig(besagBLUEs, trait = "BLUEs_yield")
## Unweighted ammi with new function (update svd in iterative process)
amNw <- gxeAmmi(besagBLUEs, trait = "BLUEs_yield")

## Identical results - no weighting and no missing values in data, so as expected.
amOrig$anova
amNw$anova

## Same as above, but number of principal components determined by algorithm.
amOrig2 <- gxeAmmiOrig(besagBLUEs, trait = "BLUEs_yield", nPC = NULL)
amNw2 <- gxeAmmi(besagBLUEs, trait = "BLUEs_yield", nPC = NULL)

## 2 PC in original function (forward selection),
## 1 PC in new function (method from Malik 2018)
amOrig2$anova
amNw2$anova

## Now using weights, 2PC
amOrigWt <- gxeAmmiOrig(besagBLUEs, trait = "BLUEs_yield", useWt = TRUE)
amNwWt <- gxeAmmi(besagBLUEs, trait = "BLUEs_yield", useWt = TRUE)

## Here results are different, but in the original function weighting means
## computing weighted residuals and the performing svd, which is very different
## from the new method.
amOrigWt$anova
amNwWt$anova

## Similar but with 5PC
amOrigWt2 <- gxeAmmiOrig(besagBLUEs, trait = "BLUEs_yield", useWt = TRUE, nPC = 5)
amNwWt2 <- gxeAmmi(besagBLUEs, trait = "BLUEs_yield", useWt = TRUE, nPC = 5)

## Problem is very obvious here. SS Residuals = -2022 for new function.
## Also SS PC1 different then SS PC1 when only 2 PC are computed.
amOrigWt2$anova
amNwWt2$anova


