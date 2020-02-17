library(statgenSTA)
library(agridat)

## Simple example used in article.
TDBesag <- createTD(besag.met, genotype = "gen", trial = "county",
                    rowCoord = "row", colCoord = "col", repId = "rep")
STABesag <- fitTD(TDBesag, design = "res.rowcol", traits = "yield")
BLUEsBesag <- STAtoTD(STABesag, what = c("BLUEs", "seBLUEs"), addWt = TRUE)

geAmBesagOrig <- gxeAmmiOrig(BLUEsBesag, trait = "BLUEs_yield")
geAmBesagOrigWt <- gxeAmmiOrig(BLUEsBesag, trait = "BLUEs_yield", useWt = TRUE)
geAmBesagNw <- gxeAmmi(BLUEsBesag, trait = "BLUEs_yield")
geAmBesagNwWt <- gxeAmmi(BLUEsBesag, trait = "BLUEs_yield", useWt = TRUE)

geGGEBesagOrig <- gxeGGEOrig(BLUEsBesag, trait = "BLUEs_yield")
geGGEBesagNw <- gxeGGE(BLUEsBesag, trait = "BLUEs_yield")
geGGEBesagNwWt <- gxeGGE(BLUEsBesag, trait = "BLUEs_yield", useWt = TRUE)

abs(geAmBesagOrig$envScores) - abs(geAmBesagNw$envscores)
abs(geGGEBesagOrig$envScores) - abs(geGGEBesagNw$envscores)

geAmBesagOrig$genoScores / geAmBesagNw$genoscores
geGGEBesagOrig$genoScores / geGGEBesagNw$genoscores

## Add some random missings to BLUEs
BLUEsBesag2 <- BLUEsBesag
BLUEsBesag2$C1$BLUEs_yield[1] <- NA
geAmBesagNw2 <- gxeAmmi(BLUEsBesag2, trait = "BLUEs_yield", useWt = TRUE)
geAmBesagNw3 <- gxeAmmi(BLUEsBesag2, trait = "BLUEs_yield", useWt = TRUE)
geAmBesagNw2$envScores - geAmBesagNw$envScores
geAmBesagNw2$genoScores - geAmBesagNw$genoScores
geAmBesagNw$fitted - geAmBesagNw2$fitted

## Larger example without missing and without wt.
geAmMaizeOrig <-  gxeAmmiOrig(TDMaize, trait = "yld")
geAmMaizeNw <-  gxeAmmi(TDMaize, trait = "yld", nPC = 6)

## Example with a lot of missing values.

datAus <- read.csv("C:/Projects/R other projects/Australia_Wheat/wheat_Australia_2008_2009_7loc.csv")

TDAus <- createTD(datAus, genotype = "geno", trial = "Trial", rowCoord = "row",
                  colCoord = "col")
STAAus <- fitTD(TDAus, design = "rowcol", traits = "yield")
BLUEsAus <- STAtoTD(STAAus, what = c("BLUEs", "seBLUEs"), addWt = TRUE)

geAmAusOrig <- gxeAmmiOrig(BLUEsAus, trait = "BLUEs_yield")
geAmAusOrigWt <- gxeAmmiOrig(BLUEsAus, trait = "BLUEs_yield", useWt = TRUE)
geAmAusNw <- gxeAmmi(BLUEsAus, trait = "BLUEs_yield")
geAmAusNwWt <- gxeAmmi(BLUEsAus, trait = "BLUEs_yield", useWt = TRUE)

geGGEAusOrig <- gxeGGEOrig(BLUEsAus, trait = "BLUEs_yield")
geGGEAusNw <- gxeGGE(BLUEsAus, trait = "BLUEs_yield")
geGGEAusNwWt <- gxeGGE(BLUEsAus, trait = "BLUEs_yield", useWt = TRUE)


dropsTD <- statgenSTA::createTD(data = dropsPheno, genotype = "Variety_ID",
                                trial = "Experiment")
testDrops <- do.call(rbind, dropsTD)$grain.yield
dim(testDrops) <- c(246, 10)
statgenGxE:::testPPB(testDrops)

dropsAmOrig <- gxeAmmiOrig(TD = dropsTD, trait = "grain.yield", nPC = NULL)
dropsAmNw <- gxeAmmi(TD = dropsTD, trait = "grain.yield", nPC = 8)




