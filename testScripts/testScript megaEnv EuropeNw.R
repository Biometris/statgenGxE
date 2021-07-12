library(statgenSTA)
library(statgenGxE)

## Read data.
datMET <- read.csv("testScripts/METsEurope.csv", na.strings = "*")

## Filter by country.
datGer <- datMET[datMET$country. == "D" &
                   datMET$trial. == "WP2" &
                   !datMET$loc. %in% c("WORRSTADT","MAGDEBURG"), ]
datDk <- datMET[datMET$country. == "DK" &
                  datMET$trial. == "NATIONAL_LIST" &
                  !datMET$loc. %in%  c("SOLLE_STED", "STORE_MERLOSE"), ]
datNl <- datMET[datMET$country. == "NL" &
                  datMET$trial. == "IP" &
                  !datMET$loc. %in% c("LELYSTAD_LIMAGRAIN", "RILLAND",
                                      "WAGENINGEN", "TOLLEBEEK", "VEENHUIZEN"), ]
datUk <- datMET[datMET$country. == "UK" &
                  datMET$trial. == "RL_TRIAL" &
                  !datMET$loc. %in% c("INCE_BLUNDELL", "WATERBEACH","CAYTHORPE",
                                      "LITTLE_STAUGHTON", "RAGNALL", "RHYND",
                                      "HOLBEACH", "WELBOURN", "LOCKERBIE",
                                      "CROFT_ON_TEES", "ELGIN","GIRTON",
                                      "BAUMBER"), ]

## Get trials per country.
trGer <- unique(datGer$env.)
trDk <- unique(datDk$env.)
trNl <- unique(datNl$env.)
trUk <- unique(datUk$env.)

## Bind together.
datTot <- rbind(datGer, datDk, datNl, datUk)

## Create TD object.
TDMet <- createTD(datTot, genotype = "geno.", trial = "env.", year = "year.",
                  loc = "loc.")

## Compute megaEnv Ger.
megaEnvGer <- gxeMegaEnvNw(TDMet, trials = trGer, trait = "yield",
                           useWinGeno = FALSE, byYear = TRUE, cutOff = 0.8)
summary(megaEnvGer)
attr(megaEnvGer$summTab, "CRDR")

## Compute megaEnv Dk.
ammiDk <- gxeAmmi(TDMet, trials = trDk, trait = "yield",
                  byYear = TRUE, nPC = NULL)

megaEnvDk <- gxeMegaEnvNw(TDMet, trials = trDk, trait = "yield",
                          useWinGeno = FALSE, byYear = TRUE, cutOff = 0.8)
summary(megaEnvDk)
attr(megaEnvDk$summTab, "CRDR")

## Compute megaEnv NL.
ammiNl <- gxeAmmi(TDMet, trials = trNl, trait = "yield",
                  byYear = TRUE, nPC = NULL)

megaEnvNl <- gxeMegaEnvNw(TDMet, trials = trNl, trait = "yield",
                          useWinGeno = FALSE, byYear = TRUE, cutOff = 0.8)
summary(megaEnvNl)
attr(megaEnvNl$summTab, "CRDR")

## Compute megaEnv Uk.
megaEnvUk <- gxeMegaEnvNw(TDMet, trials = trUk, trait = "yield",
                          useWinGeno = FALSE, byYear = TRUE, cutOff = 0.8)
summary(megaEnvUk)
attr(megaEnvUk$summTab, "CRDR")

