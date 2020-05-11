library(statenGxE)
library(agridat)

## cornelius.maize - 20 trials
datCM <- cornelius.maize
TDCM <- createTD(datCM, genotype = "gen", trial = "env")

vcCMlme <- gxeVarComp(TDCM, trait = "yield")
diagnostics(vcCMlme)
summary(vcCMlme)
vcCMlme$fullRandVC
vc(vcCMlme)
herit(vcCMlme)

vcCMasr <- gxeVarComp(TDCM, trait = "yield", engine = "asreml")
diagnostics(vcCMasr)
summary(vcCMasr)
vcCMasr$fullRandVC
vc(vcCMasr)
herit(vcCMasr)

## dasilva.maize - 9 trial, 3 replicates.
datDM <- dasilva.maize
TDDM <- createTD(datDM, genotype = "gen", trial = "env")

vcDMlme <- gxeVarComp(TDDM, trait = "yield")
diagnostics(vcDMlme)
summary(vcDMlme)
vcDMlme$fullRandVC
vc(vcDMlme)
herit(vcDMlme)

vcDMasr <- gxeVarComp(TDDM, trait = "yield", engine = "asreml")
diagnostics(vcDMasr)
summary(vcDMasr)
vcDMasr$fullRandVC
vc(vcDMasr)
herit(vcDMasr)

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
vcASlme$fullRandVC
vc(vcASlme)
herit(vcASlme)

vcASasr <- gxeVarComp(TDAS, trait = "yield", locationYear = TRUE, engine = "asreml")
diagnostics(vcASasr)
summary(vcASasr)
vcASasr$fullRandVC
vc(vcASasr)
herit(vcASasr)

## australia.soybean - 4 locations, 2 years - crossed
table(australia.soybean$year, australia.soybean$loc)

datAS2 <- australia.soybean
## Create TD object.
TDAS2 <- createTD(datAS2, genotype = "gen", trial = "env", loc = "loc",
                  year = "year")

vcAS2lme <- gxeVarComp(TDAS2, trait = "yield", locationYear = TRUE)
diagnostics(vcAS2lme)
summary(vcAS2lme)
vcAS2lme$fullRandVC
vc(vcAS2lme)
herit(vcAS2lme)

vcAS2asr <- gxeVarComp(TDAS2, trait = "yield", locationYear = TRUE, engine = "asreml")
diagnostics(vcAS2asr)
summary(vcAS2asr)
vcAS2asr$fullRandVC
vc(vcAS2asr)
herit(vcAS2asr)

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
vcSRlme$fullRandVC
vc(vcSRlme)
herit(vcSRlme)

vcSRasr <- gxeVarComp(TDSR, trait = "yield", locationYear = TRUE, engine = "asreml")
diagnostics(vcSRasr)
summary(vcSRasr)
vcSRasr$fullRandVC
vc(vcSRasr)
herit(vcSRasr)

## blackman.wheat - 7 trials, 2 treatment levels.
datBW <- blackman.wheat
TDBW <- createTD(datBW, genotype = "gen", trial = "loc")

vcBWlme <- gxeVarComp(TDBW, trait = "yield", nesting = "nitro")
diagnostics(vcBWlme)
summary(vcBWlme)
vcBWlme$fullRandVC
vc(vcBWlme)
herit(vcBWlme)

vcBWasr <- gxeVarComp(TDBW, trait = "yield", nesting = "nitro", engine = "asreml")
diagnostics(vcBWasr)
summary(vcBWasr)
vcBWasr$fullRandVC
vc(vcBWasr)
herit(vcBWasr)




