library(statgenSTA)
library(statgenGxE)

## Read data.
dat <- read.csv("W:/PROJECTS/BiomDatasetsArchive/statgenPipeline/Sunflower Argentina (internal use only)/Net_Argentina_7_years_Tesis_Abelardo.csv",
                stringsAsFactors = TRUE)
## Convert Rep to factor.
dat$Rep <- factor(dat$Rep)

## MegaEnv for original data.

dat$region <- dat$Region
dat$Env <- interaction(dat$Loc, dat$Year, sep = "_")
TDOrig <- createTD(data = dat, genotype = "Hybrid", trial = "Env",
                   loc = "Loc", year = "Year", subBlock = "Rep")
fitModOrig <- fitTD(TDOrig, what = "fixed", design = "ibd", traits = "oilyield")
BLUEsOrig <- STAtoTD(fitModOrig, what = "BLUEs",
                     keep = c("year", "loc", "region"))
megaEnvOrig <- gxeMegaEnvNw(BLUEsOrig, trait = "oilyield", useWinGeno = FALSE,
                            byYear = TRUE, cutOff = 0.8)
summary(megaEnvOrig)
attr(megaEnvOrig$summTab, "CRDR")
fullDatOrig <- do.call(rbind, megaEnvOrig$TD)
table(fullDatOrig$region, fullDatOrig$megaEnv)




## Fit variance components model.
fitMod <- asreml::asreml(fixed = oilyield ~ 1,
                         random = ~ Year + Loc + Year:Loc + Region +
                           Region:Year + Year:Loc:Rep +
                           Hybrid + Hybrid:Region +
                           Hybrid:Loc + Hybrid:Year + Hybrid:Region:Year +
                           Hybrid:Loc:Year,
                         data = dat)
summary(fitMod)

res <- list()

nYears <- 10
nLoc <- 40
nGeno <- 50

for (nReg in 2:5) {
  for (regFactor in c(1, 1.5, 2)) {
    ## Get variance components.
    vc <- summary(fitMod)$varcomp
    ## Multiply vc for genotype:region by regFactor.
    vc["Hybrid:Region", "component"] <- regFactor * vc["Hybrid:Region", "component"]

    ## Create full data.frame for simulated data.
    set.seed(1234)

    locReg <- data.frame(Loc = paste0("loc", 1:nLoc),
                         Region = sort(rep_len(1:nReg, length.out = nLoc)))

    datNew <- expand.grid(Year = paste0("year", 1:nYears),
                          Loc = paste0("loc", 1:nLoc),
                          Hybrid = paste0("hyb", 1:nGeno))
    datNew <- merge(datNew, locReg)

    noiseY <- data.frame(Year = unique(datNew$Year),
                         noiseY = rnorm(nYears, 0, vc["Year", "component"]))
    noiseL <- data.frame(Loc = unique(datNew$Loc),
                         noiseL = rnorm(nLoc, 0, vc["Loc", "component"]))
    noiseYL <- expand.grid(Year = unique(datNew$Year), Loc = unique(datNew$Loc))
    noiseYL$noiseYL <- rnorm(nYears * nLoc, 0, vc["Year:Loc", "component"])
    noiseR <- data.frame(Region = unique(datNew$Region),
                         noiseR = rnorm(nReg, 0, vc["Region", "component"]))
    noiseRY <- expand.grid(Year = unique(datNew$Year), Region = unique(datNew$Region))
    noiseRY$noiseRY <- rnorm(nReg * nYears, 0, vc["Region:Year", "component"])
    noiseG <- data.frame(Hybrid = unique(datNew$Hybrid),
                         noiseG = rnorm(nGeno, 0, vc["Hybrid", "component"]))
    noiseGR <- expand.grid(Hybrid = unique(datNew$Hybrid), Region = unique(datNew$Region))
    noiseGR$noiseGR <- rnorm(nGeno * nReg, 0, vc["Hybrid:Region", "component"])
    noiseGL <- expand.grid(Hybrid = unique(datNew$Hybrid), Loc = unique(datNew$Loc))
    noiseGL$noiseGL <- rnorm(nGeno * nLoc, 0, vc["Hybrid:Loc", "component"])
    noiseGY <- expand.grid(Hybrid = unique(datNew$Hybrid), Year = unique(datNew$Year))
    noiseGY$noiseGY <- rnorm(nGeno * nYears, 0, vc["Hybrid:Year", "component"])
    noiseGRY <- expand.grid(Hybrid = unique(datNew$Hybrid), Region = unique(datNew$Region), Year = unique(datNew$Year))
    noiseGRY$noiseGRY <- rnorm(nGeno * nReg * nYears, 0, vc["Hybrid:Region:Year", "component"])
    noiseGLY <- expand.grid(Hybrid = unique(datNew$Hybrid), Loc = unique(datNew$Loc), Year = unique(datNew$Year))
    noiseGLY$noiseGLY <- rnorm(nGeno * nLoc * nYears, 0, vc["Hybrid:Loc:Year", "component"])
    noiseRes <- rnorm(nGeno * nLoc * nYears, 0, vc["units!R", "component"])

    datNew <- Reduce(merge,
                     x = list(noiseY, noiseL, noiseYL, noiseR, noiseRY, noiseG,
                              noiseGR, noiseGL, noiseGY, noiseGRY, noiseGLY),
                     init = datNew)
    datNew$noiseRes <- noiseRes
    datNew$oilyield <- fitMod$coefficients$fixed[1, 1] +
      rowSums(datNew[, 5:ncol(datNew)])

    ## MegaEnv for new data.

    datNew$region <- datNew$Region
    datNew$Env <- interaction(datNew$Loc, datNew$Year, sep = "_")
    TDNew <- createTD(data = datNew,
                      genotype = "Hybrid", trial = "Env",
                      loc = "Loc", year = "Year")
    megaEnvNw <- gxeMegaEnvNw(TDNew, trait = "oilyield", useWinGeno = FALSE,
                              byYear = TRUE, cutOff = 0.8)
    CRDR <- attr(megaEnvNw$summTab, "CRDR")
    fullDatNw <- do.call(rbind, megaEnvNw$TD)

    res[[length(res) + 1]] <- list(dat = fullDatNw, CRDR = CRDR, nReg = nReg,
                                   regFactor = regFactor)
  }
}

sapply(res, `[[`, "nReg")
sapply(res, `[[`, "regFactor")
sapply(res, `[[`, "CRDR")
dats <- lapply(res, `[[`, "dat")
sapply(dats, FUN = function(dat) {
  length(unique(dat$megaEnv))
})
lapply(dats, FUN = function(dat) {
  table(dat$region, dat$megaEnv)
})





### Less overlap between years
## 50 hybrids in year1, 25 'old' and 25 'new' hybrids in every following year.

nYears <- 10
nLoc <- 40
nGeno <- 50 + (nYears - 1) * 25

res2 <- list()

for (nReg in 2:5) {
  for (regFactor in c(1, 1.5, 2)) {

    ## Get variance components.
    vc <- summary(fitMod)$varcomp
    ## Multiply vc for genotype:region by regFactor.
    vc["Hybrid:Region", "component"] <- regFactor * vc["Hybrid:Region", "component"]

    ## Create full data.frame for simulated data.
    set.seed(1234)

    locReg <- data.frame(Loc = paste0("loc", 1:nLoc),
                         Region = sort(rep_len(1:nReg, length.out = nLoc)))

    datNew <- expand.grid(Year = paste0("year", 1:nYears),
                          Loc = paste0("loc", 1:nLoc),
                          Hybrid = paste0("hyb", 1:nGeno))
    datNew <- merge(datNew, locReg)

    noiseY <- data.frame(Year = unique(datNew$Year),
                         noiseY = rnorm(nYears, 0, vc["Year", "component"]))
    noiseL <- data.frame(Loc = unique(datNew$Loc),
                         noiseL = rnorm(nLoc, 0, vc["Loc", "component"]))
    noiseYL <- expand.grid(Year = unique(datNew$Year), Loc = unique(datNew$Loc))
    noiseYL$noiseYL <- rnorm(nYears * nLoc, 0, vc["Year:Loc", "component"])
    noiseR <- data.frame(Region = unique(datNew$Region),
                         noiseR = rnorm(nReg, 0, vc["Region", "component"]))
    noiseRY <- expand.grid(Year = unique(datNew$Year), Region = unique(datNew$Region))
    noiseRY$noiseRY <- rnorm(nReg * nYears, 0, vc["Region:Year", "component"])
    noiseG <- data.frame(Hybrid = unique(datNew$Hybrid),
                         noiseG = rnorm(nGeno, 0, vc["Hybrid", "component"]))
    noiseGR <- expand.grid(Hybrid = unique(datNew$Hybrid), Region = unique(datNew$Region))
    noiseGR$noiseGR <- rnorm(nGeno * nReg, 0, vc["Hybrid:Region", "component"])
    noiseGL <- expand.grid(Hybrid = unique(datNew$Hybrid), Loc = unique(datNew$Loc))
    noiseGL$noiseGL <- rnorm(nGeno * nLoc, 0, vc["Hybrid:Loc", "component"])
    noiseGY <- expand.grid(Hybrid = unique(datNew$Hybrid), Year = unique(datNew$Year))
    noiseGY$noiseGY <- rnorm(nGeno * nYears, 0, vc["Hybrid:Year", "component"])
    noiseGRY <- expand.grid(Hybrid = unique(datNew$Hybrid), Region = unique(datNew$Region), Year = unique(datNew$Year))
    noiseGRY$noiseGRY <- rnorm(nGeno * nReg * nYears, 0, vc["Hybrid:Region:Year", "component"])
    noiseGLY <- expand.grid(Hybrid = unique(datNew$Hybrid), Loc = unique(datNew$Loc), Year = unique(datNew$Year))
    noiseGLY$noiseGLY <- rnorm(nGeno * nLoc * nYears, 0, vc["Hybrid:Loc:Year", "component"])
    noiseres2 <- rnorm(nGeno * nLoc * nYears, 0, vc["units!R", "component"])

    datNew <- Reduce(merge,
                     x = list(noiseY, noiseL, noiseYL, noiseR, noiseRY, noiseG,
                              noiseGR, noiseGL, noiseGY, noiseGRY, noiseGLY),
                     init = datNew)
    datNew$noiseres2 <- noiseres2
    datNew$oilyield <- fitMod$coefficients$fixed[1, 1] +
      rowSums(datNew[, 5:ncol(datNew)])

    ## Remove 25 overlapping genotypes per year.
    for (year in 1:nYears) {
      yearVal <- paste0("year", year)
      hybYear <- paste0("hyb", (25 * (year - 1) + 1):(25 * (year - 1) + 50))
      datNew <- datNew[!(datNew$Year == yearVal &
                           !datNew$Hybrid %in% hybYear), ]
    }

    ## MegaEnv for new data.

    datNew$region <- datNew$Region
    datNew$Env <- interaction(datNew$Loc, datNew$Year, sep = "_")
    TDNew <- createTD(data = datNew,
                      genotype = "Hybrid", trial = "Env",
                      loc = "Loc", year = "Year")
    megaEnvNw <- gxeMegaEnvNw(TDNew, trait = "oilyield", useWinGeno = FALSE,
                              byYear = TRUE, cutOff = 0.8)
    CRDR <- attr(megaEnvNw$summTab, "CRDR")
    fullDatNw <- do.call(rbind, megaEnvNw$TD)

    res2[[length(res2) + 1]] <- list(dat = fullDatNw, CRDR = CRDR, nReg = nReg,
                                   regFactor = regFactor)
  }
}

sapply(res2, `[[`, "nReg")
sapply(res2, `[[`, "regFactor")
sapply(res2, `[[`, "CRDR")
dats <- lapply(res2, `[[`, "dat")
sapply(dats, FUN = function(dat) {
  length(unique(dat$megaEnv))
})
lapply(dats, FUN = function(dat) {
  table(dat$region, dat$megaEnv)
})






