library(statgenSTA)
library(statgenGxE)

dat <- read.csv('C:/Projects//R other projects/Australia_Wheat/wheat_Australia_2008_2009_7loc.csv')
testTD <- createTD(dat, genotype = "geno", trial = "Trial", repId = "rep",
                   loc = "Location", year = "Year", rowCoord = "row",
                   colCoord = "col")
testMod <- fitTD(testTD, design = "res.rowcol", traits = "yield", what = "fixed")
BLUESTD <- STAtoTD(testMod, what = "BLUEs", keep = c("year", "loc"))

megaEnvOrig <- gxeMegaEnv(BLUESTD, trait = "yield")
megaEnvOrig

megaEnvNew <- gxeMegaEnvNw(BLUESTD, trait = "yield", useWinGeno = FALSE,
                           byYear = TRUE)
megaEnvNew
