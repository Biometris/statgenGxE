library(statgenSTA)
library(statgenGxE)

SB_Yield <- read.csv(system.file("extdata", "SB_yield.csv",
                                 package = "statgenSTA"),
                     stringsAsFactors = FALSE, na.strings = c("NA", "*"))
SB_Yield$year <- as.factor(substring(SB_Yield$Env, 6))
SB_Yield <- SB_Yield[SB_Yield$Genotype != "200",]

# Create object of class TD
testTD <- createTD(data = SB_Yield, genotype = "Genotype", trial = "Env",
                   repId = "Rep", subBlock = "Subblock", rowCoord = "Row",
                   colCoord = "Column", year = "year")
testSTA <- fitTD(testTD, design = "rcbd", engine = "SpATS", what = "fixed",
                 traits = "yield")
testBLUEs <- STAtoTD(testSTA, what = c("BLUEs", "seBLUEs"), addWt = TRUE,
                     keep = "year")

## Create TD object from drops data.
dropsTD <- statgenSTA::createTD(data = dropsPheno, genotype = "Variety_ID",
                                trial = "Experiment")

library(asreml)

## Basis model, just genotype and trial
vc1 <- gxeVarComp2(TD = dropsTD, trait = "grain.yield", engine = "asreml")
wald(vc1)
summary(vc1)$varcomp

vc2 <- gxeVarComp2(TD = dropsTD, trait = "grain.yield", engine = "asreml",
                  group = "scenarioWater")
wald(vc2)
summary(vc2)$varcomp

vc3 <- gxeVarComp2(TD = testTD, trait = "yield", engine = "asreml")
wald(vc3)
summary(vc3)$varcomp

vc4 <- gxeVarComp2(TD = testTD, group = "year", trait = "yield", engine = "asreml")
wald(vc4)
summary(vc4)$varcomp

