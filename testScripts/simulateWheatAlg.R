## Read data.
dat <- read.csv("W:/PROJECTS/BiomDatasetsArchive/statgenPipeline/DurumWheatAlgeria/Durum wheat Algeria.csv",
                stringsAsFactors = TRUE)
## Convert block to factor.
dat$Block <- factor(dat$Block)

## Fit variance components model.
fitMod <- asreml::asreml(fixed = yield ~ 1,
                         random = ~ Year + Location + Year:Location + Region +
                           Region:Year + Year:Location:Block + Genotype +
                           Genotype:Region + Genotype:Location +
                           Genotype:Year + Genotype:Region:Year +
                           Genotype:Location:Year,
                         data = dat)
summary(fitMod)

## Get variance components.
vc <- summary(fitMod)$varcomp
vc["Genotype:Region", "component"] <- 2 * vc["Genotype:Region", "component"]

## Create full data.frame for 20 years.
datNew <- expand.grid(Year = paste0("year", 1:20),
                      Location = unique(dat$Location),
                      Genotype = unique(dat$Genotype))
datNew <- merge(datNew, unique(dat[c("Location", "Region")]))

set.seed(1)
noiseY <- data.frame(Year = unique(datNew$Year),
                     noiseY = rnorm(20, 0, vc["Year", "component"]))
noiseL <- data.frame(Location = unique(datNew$Location),
                     noiseL = rnorm(11, 0, vc["Location", "component"]))
noiseYL <- expand.grid(Year = unique(datNew$Year), Location = unique(datNew$Location))
noiseYL$noiseYL <- rnorm(20 * 11, 0, vc["Year:Location", "component"])
noiseR <- data.frame(Region = unique(datNew$Region),
                     noiseR = rnorm(2, 0, vc["Region", "component"]))
noiseRY <- expand.grid(Year = unique(datNew$Year), Region = unique(datNew$Region))
noiseRY$noiseRY <- rnorm(2 * 20, 0, vc["Region:Year", "component"])
noiseG <- data.frame(Genotype = unique(datNew$Genotype),
                     noiseG = rnorm(24, 0, vc["Genotype", "component"]))
noiseGR <- expand.grid(Genotype = unique(datNew$Genotype), Region = unique(datNew$Region))
noiseGR$noiseGR <- rnorm(24 * 2, 0, vc["Genotype:Region", "component"])
noiseGL <- expand.grid(Genotype = unique(datNew$Genotype), Location = unique(datNew$Location))
noiseGL$noiseGL <- rnorm(24 * 11, 0, vc["Genotype:Location", "component"])
noiseGY <- expand.grid(Genotype = unique(datNew$Genotype), Year = unique(datNew$Year))
noiseGY$noiseGY <- rnorm(24 * 20, 0, vc["Genotype:Year", "component"])
noiseGRY <- expand.grid(Genotype = unique(datNew$Genotype), Region = unique(datNew$Region), Year = unique(datNew$Year))
noiseGRY$noiseGRY <- rnorm(24 * 2 * 20, 0, vc["Genotype:Region:Year", "component"])
noiseGLY <- expand.grid(Genotype = unique(datNew$Genotype), Location = unique(datNew$Location), Year = unique(datNew$Year))
noiseGLY$noiseGLY <- rnorm(24 * 11 * 20, 0, vc["Genotype:Location:Year", "component"])
noiseRes <- rnorm(24 * 11 * 20, 0, vc["units!R", "component"])

datNew <- Reduce(merge,
                 x = list(noiseY, noiseL, noiseYL, noiseR, noiseRY, noiseG,
                          noiseGR, noiseGL, noiseGY, noiseGRY, noiseGLY),
                 init = datNew)
datNew$noiseRes <- noiseRes
datNew$yield <- fitMod$coefficients$fixed[1, 1] + rowSums(datNew[, 5:ncol(datNew)])

## MegaEnv for original data.

dat$region <- dat$Region
TDOrig <- createTD(data = dat, genotype = "Genotype", trial = "Environment",
                   loc = "Location", year = "Year", subBlock = "Block")
fitModOrig <- fitTD(TDOrig, what = "fixed", design = "ibd", traits = "yield")
BLUEsOrig <- STAtoTD(fitModOrig, what = "BLUEs",
                     keep = c("year", "loc", "region"))
megaEnvOrig <- gxeMegaEnvNw(BLUEsOrig, trait = "yield", useWinGeno = FALSE,
                            byYear = TRUE, cutOff = 0.8)
summary(megaEnvOrig)
fullDatOrig <- do.call(rbind, megaEnvOrig$TD)
table(fullDatOrig$region, fullDatOrig$megaEnv)

## MegaEnv for new data.

datNew$region <- datNew$Region
datNew$Environment <- interaction(datNew$Location, datNew$Year, sep = "_")
TDNew <- createTD(data = datNew,
                  genotype = "Genotype", trial = "Environment",
                   loc = "Location", year = "Year")
megaEnvNw <- gxeMegaEnvNw(TDNew, trait = "yield", useWinGeno = FALSE,
                          byYear = TRUE, cutOff = 0.8)
summary(megaEnvNw)
attr(megaEnvNw$summTab, "CRDR")


