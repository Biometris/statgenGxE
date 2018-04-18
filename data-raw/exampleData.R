## Create TDHeat05.

# Read raw data
SB_Yield <- read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
                     stringsAsFactors = FALSE, na.strings = c("NA", "*"))
# Restrict to HEAT05
Heat05 <- SB_Yield[SB_Yield$Env == "HEAT05", ]
# Create object of class TD
TDHeat05 <- createTD(data = Heat05, genotype = "Genotype", trial = "Env",
                     repId = "Rep", subBlock = "Subblock", rowId = "Row",
                     colId = "Column", rowCoord = "Row",
                     colCoord = "Column")
# Export to package
devtools::use_data(TDHeat05, overwrite = TRUE)

## Create TDMaize.

# Read raw data
F2Maize <- read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
                    stringsAsFactors = FALSE)
# Create object of class TD
TDMaize <- createTD(data = F2Maize, genotype = "genotype.", trial = "env.")
# Export to package
devtools::use_data(TDMaize, overwrite = TRUE)

## Create a dataset for unit testing.
set.seed(123)
testData <- data.frame(seed = rep(x = paste0("G", rep(x = 1:15, times = 2)), times = 3),
                       family = paste0("F", 1:90),
                       field = rep(x = paste0("E", 1:3), each = 30),
                       rep = rep(x = c(1, 2), each = 15),
                       checkId = sample.int(n = 2, size = 90, replace = TRUE),
                       X = rep(x = rep(x = 1:3, each = 10), times = 3),
                       Y = as.numeric(replicate(n = 9, expr = sample.int(n = 10))),
                       block = rep(x = 1:5, times = 18),
                       t1 = rnorm(n = 90, mean = 80, sd = 25),
                       t2 = rnorm(n = 90, mean = 3, sd = 0.5),
                       t3 = rnorm(n = 90, mean = 60, sd = 10),
                       t4 = rnorm(n = 90, mean = 80, sd = 25)
)
## Add some random NAs to traits t3 and t4.
testData$t3[sample.int(n = 90, size = 15)] <- NA
testData$t4[sample.int(n = 90, size = 15)] <- NA
## Export to package
devtools::use_data(testData, overwrite = TRUE)

## Use data from qtl package for testing cross functions.
## Save locally to prevent errors from changes in data.
data(fake.f2, package = "qtl")
## Restrict data.
testF2 <- fake.f2[c(1,2,"X"), 1:50]
## Add a duplicate marker on chr 1 for testing purposes.
testF2$geno$`1`$data <- cbind(testF2$geno$`1`$data, testF2$geno$`1`$data[, 7])
colnames(testF2$geno$`1`$data)[8] <- "D1M37"
testF2$geno$`1`$map <- c(testF2$geno$`1`$map, 137.37)
names(testF2$geno$`1`$map)[8] <- "D1M37"
## Export to package
devtools::use_data(testF2, overwrite = TRUE)

data(fake.4way, package = "qtl")
## Restrict data.
test4way <- fake.4way[c(1,2,"X"), 1:50]
## Export to package
devtools::use_data(test4way, overwrite = TRUE)

data(fake.bc, package = "qtl")
## Restrict data.
testBc <- fake.bc[1:3, 1:50]
## Export to package.
devtools::use_data(testBc, overwrite = TRUE)

## Create data for vignette.
# Read raw data.
wheatAus <- read.csv(system.file("extdata", "wheat_Australia_recoded.csv",
                                 package = "RAP"), stringsAsFactors = FALSE)
wheatAus$year <- 2000 + as.numeric(substring(text = wheatAus$Trial, first = 1,
                                             last = 2))
wheatAus$loc <- substring(text = wheatAus$Trial, first = 3)
# Export to package.
devtools::use_data(wheatAus, overwrite = TRUE)

# Read raw data.
dat2011 <- read.delim(system.file("extdata", "pheno_data2011.txt",
                                  package = "RAP"))
dat2012 <- read.delim(system.file("extdata", "pheno_data2012.txt",
                                  package = "RAP"))
# Split data into separate year/trials.
dat2011_1 <- dat2011[, 1:11]
dat2011_1$trial <- "SR_FI_11"
colnames(dat2011_1)[8:11] <- c("DH", "GY", "NKS", "TKW")
dat2011_2 <- dat2011[, c(1:7, 12:15)]
dat2011_2$trial <- "SR_MWS_11"
colnames(dat2011_2)[8:11] <- c("DH", "GY", "NKS", "TKW")
dat2011tot <- rbind(dat2011_1, dat2011_2)
dat2011tot$year <- 2011

dat2012_1 <- dat2012[, 1:8]
dat2012_1$trial <- "SR_FI_12"
colnames(dat2012_1)[8] <- "GY"
dat2012_2 <- dat2012[, c(1:7, 9)]
dat2012_2$trial <- "SR_MWS_12"
colnames(dat2012_2)[8] <- "GY"
dat2012_3 <- dat2012[, c(1:7, 10)]
dat2012_3$trial <- "C_SWS_12"
colnames(dat2012_3)[8] <- "GY"
dat2012tot <- rbind(dat2012_1, dat2012_2, dat2012_3)
dat2012tot$year <- 2012
dat2012tot[c("DH", "NKS", "TKW")] <- NA
# Bind year data together and rename genotypes.
wheatChl <- rbind(dat2011tot, dat2012tot)
wheatChl$trt_id <- sprintf("G%03d", wheatChl$trt_id)
# Export to package
devtools::use_data(wheatChl, overwrite = TRUE)







