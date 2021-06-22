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
                           useWinGeno = FALSE, byYear = TRUE)
summary(megaEnvGer)
attr(megaEnvGer$summTab, "CRDR")

## Compute megaEnv Dk.
ammiDk <- gxeAmmi(TDMet, trials = trDk, trait = "yield",
                  byYear = TRUE, nPC = NULL)

megaEnvDk <- gxeMegaEnvNw(TDMet, trials = trDk, trait = "yield",
                          useWinGeno = FALSE, byYear = TRUE)
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
                          useWinGeno = FALSE, byYear = TRUE)
summary(megaEnvUk)
attr(megaEnvUk$summTab, "CRDR")

