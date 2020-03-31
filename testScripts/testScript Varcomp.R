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
vc1 <- gxeVarComp(TD = dropsTD1, trait = "grain.yield", engine = "asreml")
vc1$call$fixed
vc1$call$random
wald(vc1)
summary(vc1)$varcomp

## Add a group variable to the model.
vc2 <- gxeVarComp(TD = dropsTD1, trait = "grain.yield", engine = "asreml",
                   trialGroup = "scenarioWater")
vc2$call$fixed
vc2$call$random
wald(vc2)
summary(vc2)$varcomp

## Basic model, just genotype and loc x year - drops data.
vc1a <- gxeVarComp(TD = dropsTD2, trait = "grain.yield", engine = "asreml")
vc1a$call$fixed
vc1a$call$random
wald(vc1a)
summary(vc1a)$varcomp

## Add a group variable to the model - for loc x year.
vc2a <- gxeVarComp(TD = dropsTD2, trait = "grain.yield", engine = "asreml",
                    trialGroup = "scenarioWater")
vc2a$call$fixed
vc2a$call$random
wald(vc2a)
summary(vc2a)$varcomp

## Basic model, just genotype and loc x year - complete data set.
vc1b <- gxeVarComp(TD = dropsTD3, trait = "grain.yield", engine = "asreml")
vc1b$call$fixed
vc1b$call$random
wald(vc1b)
summary(vc1b)$varcomp

## Add a group variable to the model - for loc x year.
vc2b <- gxeVarComp(TD = dropsTD3, trait = "grain.yield", engine = "asreml",
                    trialGroup = "scenarioWater")
vc2b$call$fixed
vc2b$call$random
wald(vc2b)
summary(vc2b)$varcomp

## Basic model for replicated data.
vc3 <- gxeVarComp(TD = testTD, trait = "yield", engine = "asreml")
vc3$call$fixed
vc3$call$random
wald(vc3)
summary(vc3)$varcomp

## Add a group variable - using year as group here.
vc4 <- gxeVarComp(TD = testTD, trialGroup = "year", trait = "yield",
                   engine = "asreml")
vc4$call$fixed
vc4$call$random
wald(vc4)
summary(vc4)$varcomp

