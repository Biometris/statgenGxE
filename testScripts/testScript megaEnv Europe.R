library(statgenSTA)
library(statgenGxE)

## Read data.
datGer <- read.csv("w:/PROJECTS/BiomDatasetsArchive/statgenPipeline/MET_Europe/WheatYieldD.csv",
                   quote = "'", na.strings = "*")
datDk <- read.csv("w:/PROJECTS/BiomDatasetsArchive/statgenPipeline/MET_Europe/WheatYieldDK.csv",
                   quote = "'", na.strings = "*")
datNl <- read.csv("w:/PROJECTS/BiomDatasetsArchive/statgenPipeline/MET_Europe/WheatYieldNL.csv",
                   quote = "'", na.strings = "*")
datUk <- read.csv("w:/PROJECTS/BiomDatasetsArchive/statgenPipeline/MET_Europe/WheatYieldUK.csv",
                   quote = "'", na.strings = "*")

## Filter for locations used in paper.
# Names don't match at all.
datGer <- datGer[datGer$loc %in%
                   c("DACHWIG", "FUTTERKAMP_1", "GIEBELSTADT", "HASSLOCH",
                     "KERPEN-BUIR_1", "KONIGSLUTTER_1", "LADENBURG", "NEUHOF_1",
                     "NEUHOF_2", "NOSSEN", "OLVENSTEDT", "OSTINGHAUSEN",
                     "PRENZLAU", "RAUISCHHOLZHAUSEN", "RETHMAR",
                     "SOPHIENHOF_1", "VEINAU_4"), ]
# SOLL = HINNE_RUP ???
datDk <- datDk[datDk$loc %in%
                 c("GRIND_STED", "HOLEBY", "HOLSTE_BRO", "HORSENS",
                   "NYKOBING_F", "ORBAEK", "SKAELSKOR", "SONDER_BORG",
                   "TOLLOSE", "TRIGE", "VRA", "HINNE_RUP"), ]
datNl <- datNl[datNl$loc %in%
                 c("DRONTEN", "ZELDER", "COLIJNSPLAAT", "HORNHUIZEN",
                   "LELYSTAD_(PPO)", "NIEUW_BEERTA", "SIBCULO", "WESTMAAS",
                   "WIJNANDSRADE", "ANGEREN", "VEENHUIZEN") &
                 datNl$year > 2005, ]
# Names don't match at all.
datUk <- datUk[datUk$loc %in%
                 c("BARNSTON", "BEAUSALE", "BOWSDEN", "CAMBRIDGE", "DIDBROOK",
                   "FRAMLINGHAM", "FRISBY", "HAMPSHIRE_UP_SOMBORNE", "HASELEY",
                   "HISTON", "HUMBIE", "HYTHE", "LAURENCEKIRK", "LIMAVADY", "LYMN",
                   "PERTH", "RUDSTON", "SPALDING", "TELFORD", "ULCEBY", "WARWICK",
                   "WEST_CHARLETON", "WINCHESTER", "WOLFERTON"), ]

## Get trials per country.
trGer <- unique(datGer$env)
trDk <- unique(datDk$env)
trNl <- unique(datNl$env)
trUk <- unique(datUk$env)

## Bind together.
datTot <- rbind(datGer, datDk, datNl, datUk)

## Create TD object.
TDMet <- createTD(datTot, genotype = "geno", trial = "env")

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

