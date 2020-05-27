library(statgenSTA)
library(statgenGxE)
library(agridat)

## cornelius.maize - 20 trials
datCM <- cornelius.maize
TDCM <- createTD(datCM, genotype = "gen", trial = "env")

vcCMlme <- gxeVarComp(TDCM, trait = "yield")
diagnostics(vcCMlme)
summary(vcCMlme)
plot(vcCMlme)
vc(vcCMlme)
herit(vcCMlme)
predict(vcCMlme)

vcCMasr <- gxeVarComp(TDCM, trait = "yield", engine = "asreml")
diagnostics(vcCMasr)
summary(vcCMasr)
plot(vcCMasr)
vc(vcCMasr)
herit(vcCMasr)
predict(vcCMasr)

## dasilva.maize - 9 trial, 3 replicates.
datDM <- dasilva.maize
TDDM <- createTD(datDM, genotype = "gen", trial = "env")

vcDMlme <- gxeVarComp(TDDM, trait = "yield")
diagnostics(vcDMlme)
summary(vcDMlme)
plot(vcDMlme)
vc(vcDMlme)
herit(vcDMlme)
predict(vcDMlme)

vcDMasr <- gxeVarComp(TDDM, trait = "yield", engine = "asreml")
diagnostics(vcDMasr)
summary(vcDMasr)
plot(vcDMasr)
vc(vcDMasr)
herit(vcDMasr)
predict(vcDMasr)

## adugna.sorghum - 3 locations, 5 years - nested
table(adugna.sorghum$year, adugna.sorghum$loc)

datAS <- adugna.sorghum
## Rename existing trial column.
datAS$trialOrig <- datAS$trial
datAS$trial <- NULL
## Create TD object.
TDAS <- createTD(datAS, genotype = "gen", trial = "env", loc = "loc",
                 year = "year")

vcASlme <- gxeVarComp(TDAS, trait = "yield", locationYear = TRUE)
diagnostics(vcASlme)
summary(vcASlme)
plot(vcASlme)
vc(vcASlme)
herit(vcASlme)
predict(vcASlme)
predict(vcASlme, predictLevel = "trial")

vcASasr <- gxeVarComp(TDAS, trait = "yield", locationYear = TRUE, engine = "asreml")
diagnostics(vcASasr)
summary(vcASasr)
plot(vcASasr)
vc(vcASasr)
herit(vcASasr)
predict(vcASasr)

## australia.soybean - 4 locations, 2 years - crossed
table(australia.soybean$year, australia.soybean$loc)

datAS2 <- australia.soybean
## Create TD object.
TDAS2 <- createTD(datAS2, genotype = "gen", trial = "env", loc = "loc",
                  year = "year")

vcAS2lme <- gxeVarComp(TDAS2, trait = "yield", locationYear = TRUE)
diagnostics(vcAS2lme)
summary(vcAS2lme)
plot(vcAS2lme)
vc(vcAS2lme)
herit(vcAS2lme)
predict(vcAS2lme)

vcAS2asr <- gxeVarComp(TDAS2, trait = "yield", locationYear = TRUE, engine = "asreml")
diagnostics(vcAS2asr)
summary(vcAS2asr)
plot(vcAS2asr)
vc(vcAS2asr)
herit(vcAS2asr)
predict(vcAS2asr)

## shafii.rapeseed - 14 locations, 3 years, 3 replicates
datSR <- shafii.rapeseed
## Add trial to data.
datSR$trial <- interaction(datSR$loc, datSR$year)

## Create TD object.
TDSR <- createTD(datSR, genotype = "gen", trial = "trial", loc = "loc",
                 year = "year")

vcSRlme <- gxeVarComp(TDSR, trait = "yield", locationYear = TRUE)
diagnostics(vcSRlme)
summary(vcSRlme)
plot(vcSRlme)
vc(vcSRlme)
herit(vcSRlme)
predict(vcSRlme)

vcSRasr <- gxeVarComp(TDSR, trait = "yield", locationYear = TRUE, engine = "asreml")
diagnostics(vcSRasr)
summary(vcSRasr)
plot(vcSRasr)
vc(vcSRasr)
herit(vcSRasr)
predict(vcSRasr)
predict(vcSRasr, predictLevel = "trial")

## blackman.wheat - 7 trials, 2 treatment levels.
datBW <- blackman.wheat
TDBW <- createTD(datBW, genotype = "gen", trial = "loc")

vcBWlme <- gxeVarComp(TDBW, trait = "yield", nestingFactor = "nitro")
diagnostics(vcBWlme)
summary(vcBWlme)
plot(vcBWlme)
vc(vcBWlme)
herit(vcBWlme)
predict(vcBWlme)
predict(vcBWlme, predictLevel = "nitro")

vcBWasr <- gxeVarComp(TDBW, trait = "yield", nesting = "nitro", engine = "asreml")
diagnostics(vcBWasr)
summary(vcBWasr)
plot(vcBWasr)
vcBWasr$fullRandVC
vc(vcBWasr)
herit(vcBWasr)
predict(vcBWasr)
predict(vcBWasr, predictLevel = "nitro")

