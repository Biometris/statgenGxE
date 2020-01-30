## Create TDMaize.

# Read raw data
F2Maize <- read.csv(system.file("extdata", "F2maize_pheno.csv",
                                package = "statgenGxE"),
                    stringsAsFactors = FALSE)
# Create object of class TD
TDMaize <- createTD(data = F2Maize, genotype = "genotype.", trial = "env.")
# Export to package
usethis::use_data(TDMaize, overwrite = TRUE)

## Create a dataset for unit testing.
set.seed(1)
testData <- data.frame(seed = rep(x = paste0("G", rep(x = 1:15, times = 2)),
                                  times = 3),
                       family = rep(x = paste0("F", rep(x = 1:3, each = 5)),
                                   times = 6),
                       field = rep(x = paste0("E", 1:3), each = 30),
                       regime = rep(x = c("D", "D", "W"), each = 30),
                       rep = rep(x = c(1, 2), each = 15),
                       checkId = sample.int(n = 2, size = 90, replace = TRUE),
                       X = rep(x = rep(x = 1:3, each = 10), times = 3),
                       Y = as.numeric(replicate(n = 9,
                                                expr = sample.int(n = 10))),
                       block = rep(x = 1:5, times = 18),
                       t1 = rnorm(n = 90, mean = 80, sd = 25),
                       t2 = rnorm(n = 90, mean = 3, sd = 0.5),
                       t3 = rnorm(n = 90, mean = 60, sd = 10),
                       t4 = rnorm(n = 90, mean = 80, sd = 25))
## Add some random NAs to traits t3 and t4.
testData[sample.int(n = 90, size = 15), "t3"] <- NA
testData[sample.int(n = 90, size = 15), "t4"] <- NA

## Create TD object from testData.
testTD <- statgenSTA::createTD(data = testData, genotype = "seed",
                               trial = "field", rowCoord = "Y", colCoord = "X")

## Create a second dataset with year info.
## Just two copies of testData with a year column added.
testDataYear <- rbind(testData, testData)
testDataYear[["field"]] = rep(x = paste0("E", 1:6), each = 30)
testDataYear[["year"]] <- rep(c(1, 2), each = 90)
testDataYear[91:180, "t1"] <- runif(90) * mean(testDataYear[["t1"]])

## Create TD object from testDataYear.
testTDYear <- statgenSTA::createTD(data = testDataYear, genotype = "seed",
                                   trial = "field", rowCoord = "Y",
                                   colCoord = "X")

## Fit models on testTD and testTDYear and extract BLUEs
modelSp <- statgenSTA::fitTD(testTD, design = "rowcol", traits = c("t1", "t2"))
BLUEs <- statgenSTA::STAtoTD(modelSp, what = "BLUEs",
                             keep = c("family", "regime"))

modelSpYear <- statgenSTA::fitTD(testTDYear, design = "rowcol",
                                 traits = c("t1", "t2"))
BLUEsYear <- statgenSTA::STAtoTD(modelSpYear, what = "BLUEs",
                                 keep = c("family", "regime", "year"))

## Export all internal data in one go to package.
usethis::use_data(testTD, testTDYear, BLUEs, BLUEsYear, overwrite = TRUE,
                  internal = TRUE)
