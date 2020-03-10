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
geGGEBesagNwWt2 <- gxeGGE(BLUEsBesag, trait = "BLUEs_yield", useWt = TRUE, nPC = 3)

abs(geAmBesagOrig$envScores) - abs(geAmBesagNw$envScores)
abs(geGGEBesagOrig$envScores) - abs(geGGEBesagNw$envScores)

geAmBesagOrig$genoScores / geAmBesagNw$genoScores
geGGEBesagOrig$genoScores / geGGEBesagNw$genoScores

testBesag <- do.call(rbind, BLUEsBesag)$BLUEs_yield
dim(testBesag) <- c(64, 6)
statgenGxE:::testPPB(testBesag)

testBesagWt <- do.call(rbind, BLUEsBesag)$BLUEs_yield * do.call(rbind, BLUEsBesag)$wt
dim(testBesagWt) <- c(64, 6)
statgenGxE:::testPPB(testBesagWt)

## Add some random missings to BLUEs
BLUEsBesag2 <- BLUEsBesag
BLUEsBesag2$C1$BLUEs_yield[5] <- NA
geAmBesagNw2 <- gxeAmmi(BLUEsBesag2, trait = "BLUEs_yield", useWt = TRUE)
geAmBesagNw2$envScores - geAmBesagNw$envScores
geAmBesagNw2$genoScores - geAmBesagNw3$genoScores
geAmBesagNw2$fitted$fittedValue - geAmBesagNwWt$fitted$fittedValue

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
geAmAusOrig3 <- gxeAmmiOrig(BLUEsAus, trait = "BLUEs_yield", nPC = 3)
geAmAusOrigWt <- gxeAmmiOrig(BLUEsAus, trait = "BLUEs_yield", useWt = TRUE)
geAmAusNw <- gxeAmmi(BLUEsAus, trait = "BLUEs_yield") #376
geAmAusNwWt <- gxeAmmi(BLUEsAus, trait = "BLUEs_yield", useWt = TRUE) #4715

geGGEAusOrig <- gxeGGEOrig(BLUEsAus, trait = "BLUEs_yield")
geGGEAusNw <- gxeGGE(BLUEsAus, trait = "BLUEs_yield")
geGGEAusNwWt <- gxeGGE(BLUEsAus, trait = "BLUEs_yield", useWt = TRUE)


dropsTD <- statgenSTA::createTD(data = dropsPheno, genotype = "Variety_ID",
                                trial = "Experiment")
testDrops <- do.call(rbind, dropsTD)$grain.yield
dim(testDrops) <- c(246, 10)
statgenGxE:::testPPB(testDrops)

dropsAmOrig <- gxeAmmiOrig(TD = dropsTD, trait = "grain.yield", nPC = NULL)
dropsAmNw <- gxeAmmi(TD = dropsTD, trait = "grain.yield", nPC = 7)




