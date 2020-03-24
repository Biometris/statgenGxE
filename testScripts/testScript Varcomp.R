library(statgenSTA)
library(statgenGxE)
library(asreml)

## Read data from statgenSTA - 4 trials, 169 genotypes.
SB_Yield <- read.csv(system.file("extdata", "SB_yield.csv",
                                 package = "statgenSTA"),
                     stringsAsFactors = FALSE, na.strings = c("NA", "*"))
## Add a year variable.
SB_Yield$year <- as.factor(substring(SB_Yield$Env, 6))
## Genotype 200 only appears in one trial - remove it.
SB_Yield <- SB_Yield[SB_Yield$Genotype != "200",]

# Create object of class TD
testTD <- createTD(data = SB_Yield, genotype = "Genotype", trial = "Env",
                   repId = "Rep", subBlock = "Subblock", rowCoord = "Row",
                   colCoord = "Column", year = "year")
## Fit models en compute BLUEs.
testSTA <- fitTD(testTD, design = "rcbd", engine = "SpATS", what = "fixed",
                 traits = "yield")
testBLUEs <- STAtoTD(testSTA, what = c("BLUEs", "seBLUEs"), addWt = TRUE,
                     keep = "year")

## Create TD object from drops data from statgenGxE
dropsTD <- statgenSTA::createTD(data = dropsPheno, genotype = "Variety_ID",
                                trial = "Experiment")

## Basic model, just genotype and trial - drops data.
vc1 <- gxeVarComp2(TD = dropsTD, trait = "grain.yield", engine = "asreml")
vc1$call$fixed
vc1$call$random
wald(vc1)
summary(vc1)$varcomp

## Add a group variable to the model.
vc2 <- gxeVarComp2(TD = dropsTD, trait = "grain.yield", engine = "asreml",
                  group = "scenarioWater")
vc2$call$fixed
vc2$call$random
wald(vc2)
summary(vc2)$varcomp

## Basic model for replicated data.
vc3 <- gxeVarComp2(TD = testTD, trait = "yield", engine = "asreml")
vc3$call$fixed
vc3$call$random
wald(vc3)
summary(vc3)$varcomp

## Add a group variable - using year as group here.
vc4 <- gxeVarComp2(TD = testTD, group = "year", trait = "yield",
                   engine = "asreml")
vc4$call$fixed
vc4$call$random
wald(vc4)
summary(vc4)$varcomp

