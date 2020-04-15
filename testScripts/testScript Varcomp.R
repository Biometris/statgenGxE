## Create data set from drops with 4 year x loc combinations.
{
  ## Load data from zip file.
  phenoFile <- system.file("extdata", "grainYield_components_BLUEs.zip",
                           package = "statgenGxE")
  dropsPhenoTot  <- read.csv(unz(description = phenoFile,
                                 filename = "2b-GrainYield_components_BLUEs_level.csv"))

  ## Add year column.
  dropsPhenoTot [["year"]] <- paste0("20", substring(dropsPhenoTot [["Experiment"]],
                                                     first = 4, last = 5))
  dropsPhenoTot [["scenarioWater"]] <- factor(substr(dropsPhenoTot$Experiment,6,6))


  dropsPhenoTot$location = substr(dropsPhenoTot$Experiment, 1, 3)

  dropsPhenoTot <- droplevels(dropsPhenoTot)
  dropsPhenoTot <- dropsPhenoTot[dropsPhenoTot$location %in% c("Cam", "Deb", "Kar", "Ner"),]
  rm(list = "phenoFile")
}

library(statgenSTA)
library(statgenGxE)
library(asreml)
library(lme4)

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

## Create TD objects from drops data from statgenGxE.
dropsPheno$location = substr(dropsPheno$Experiment, 1, 3)
## First one without location.
dropsTD1 <- statgenSTA::createTD(data = dropsPheno, genotype = "Variety_ID",
                                 trial = "Experiment")
## Second one with location so year x location can be used.
dropsTD2 <- statgenSTA::createTD(data = dropsPheno, genotype = "Variety_ID",
                                 trial = "Experiment", loc = "location")
## Third ond based on 4 trials that are present in all 4 years.
dropsTD3 <- statgenSTA::createTD(data = dropsPhenoTot, genotype = "Variety_ID",
                                 trial = "Experiment", loc = "location")

## Basic model, just genotype and trial - drops data.
vc0 <- gxeVarComp(TD = dropsTD1, trait = "grain.yield", engine = "lme4")
vc0$fitMod@call$formula
vc(vc0)
herit(vc0)
p0 <- predict(vc0)

vc1 <- gxeVarComp(TD = dropsTD1, trait = "grain.yield", engine = "asreml")
vc1$fitMod$call$fixed
vc1$fitMod$call$random
wald(vc1$fitMod)
vc(vc1)
herit(vc1)
p1 <- predict(vc1)
p1_1 <- predict(vc1, predictLevel = "trial")

## Add a group variable to the model.
vc20 <- gxeVarComp(TD = dropsTD1, trait = "grain.yield", engine = "lme4",
                   trialGroup = "scenarioWater")

attr(getME(vc20$fitMod,"X"),"col.dropped")

vc20$fitMod@call$formula
vc(vc20)
herit(vc20)
p20 <- predict(vc20)

vc2 <- gxeVarComp(TD = dropsTD1, trait = "grain.yield", engine = "asreml",
                  trialGroup = "scenarioWater")
vc2$fitMod$call$fixed
vc2$fitMod$call$random
wald(vc2$fitMod)
vc(vc2)
herit(vc2)
p2 <- predict(vc2)
p2_1 <- predict(vc2, predictLevel = "trial")
p2_2 <- predict(vc2, predictLevel = "trialGroup")

## Basic model, just genotype and loc x year - drops data.
vc1a0 <- gxeVarComp(TD = dropsTD2, trait = "grain.yield", engine = "lme4")
vc1a0$fitMod@call$formula
vc(vc1a0)
herit(vc1a0)

vc1a <- gxeVarComp(TD = dropsTD2, trait = "grain.yield", engine = "asreml")
vc1a$fitMod$call$fixed
vc1a$fitMod$call$random
wald(vc1a$fitMod)
vc(vc1a)
herit(vc1a)
p1a <- predict(vc1a)
p1a_1 <- predict(vc1a, predictLevel = "trial")

## Add a group variable to the model - for loc x year.
vc2a0 <- gxeVarComp(TD = dropsTD2, trait = "grain.yield", engine = "lme4",
                    trialGroup = "scenarioWater")

attr(getME(vc2a0$fitMod,"X"),"col.dropped")

vc(vc2a0)
herit(vc2a0)

vc2a <- gxeVarComp(TD = dropsTD2, trait = "grain.yield", engine = "asreml",
                   trialGroup = "scenarioWater")
vc2a$fitMod$call$fixed
vc2a$fitMod$call$random
wald(vc2a$fitMod)
vc(vc2a)
herit(vc2a)
p2a <- predict(vc2a)
p2a_1 <- predict(vc2a, predictLevel = "trial")
p2a_2 <- predict(vc2a, predictLevel = "trialGroup")

## Basic model, just genotype and loc x year - complete data set.
vc1b0 <- gxeVarComp(TD = dropsTD3, trait = "grain.yield", engine = "lme4")
vc(vc1b0)
herit(vc1b0)
p1b0 <- predict(vc1b0)

vc1b <- gxeVarComp(TD = dropsTD3, trait = "grain.yield", engine = "asreml")
vc1b$fitMod$call$fixed
vc1b$fitMod$call$random
wald(vc1b$fitMod)
vc(vc1b)
herit(vc1b)
p1b <- predict(vc1b)
p1b_1 <- predict(vc1b, predictLevel = "trial")

## Add a group variable to the model - for loc x year.
vc2b0 <- gxeVarComp(TD = dropsTD3, trait = "grain.yield", engine = "lme4",
                    trialGroup = "scenarioWater")
vc2b0$fitMod@call$formula
vc(vc2b0)
herit(vc2b0)

vc2b <- gxeVarComp(TD = dropsTD3, trait = "grain.yield", engine = "asreml",
                   trialGroup = "scenarioWater")
vc2b$fitMod$call$fixed
vc2b$fitMod$call$random
wald(vc2b$fitMod)
vc(vc2b)
herit(vc2b)
p2b <- predict(vc2b)
p2b_1 <- predict(vc2b, predictLevel = "trial")
p2b_2 <- predict(vc2b, predictLevel = "trialGroup")

## Basic model for replicated data.
vc30 <- gxeVarComp(TD = testTD, trait = "yield", engine = "lme4")
vc(vc30)
herit(vc30)

vc3 <- gxeVarComp(TD = testTD, trait = "yield", engine = "asreml")
vc3$fitMod$call$fixed
vc3$fitMod$call$random
wald(vc3$fitMod)
vc(vc3)
herit(vc3)
p3 <- predict(vc3)
p3_1 <- predict(vc3, predictLevel = "trial")

## Add a group variable - using year as group here.
vc40 <- gxeVarComp(TD = testTD, trialGroup = "year", trait = "yield",
                   engine = "lme4")
vc(vc40)
herit(vc40)

vc4 <- gxeVarComp(TD = testTD, trialGroup = "year", trait = "yield",
                  engine = "asreml")
vc4$fitMod$call$fixed
vc4$fitMod$call$random
wald(vc4$fitMod)
vc(vc4)
herit(vc4)
p4 <- predict(vc4)
p4_1 <- predict(vc4, predictLevel = "trial")
p4_2 <- predict(vc4, predictLevel = "trialGroup")

