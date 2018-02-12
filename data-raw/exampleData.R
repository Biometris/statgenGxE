## Create TDHeat05.

# Read raw data
SB_Yield <- read.csv(system.file("extdata", "SB_yield.csv", package = "RAP"),
                     stringsAsFactors = FALSE, na.strings = c("NA", "*"))
# Restrict to HEAT05
Heat05 <- SB_Yield[SB_Yield$Env == "HEAT05", ]
# Create object of class TD
TDHeat05 <- createTD(data = Heat05, genotype = "Genotype", env = "Env",
                     repId = "Rep", subBlock = "Subblock", rowId = "Row",
                     colId = "Column", rowCoordinates = "Row",
                     colCoordinates = "Column")
# Export to package
devtools::use_data(TDHeat05, overwrite = TRUE)

## Create TDMaize

# Read raw data
F2Maize <- read.csv(system.file("extdata", "F2maize_pheno.csv", package = "RAP"),
                    stringsAsFactors = FALSE)
# Create object of class TD
TDMaize <- createTD(data = F2Maize, genotype = "genotype.", env = "env.")
# Export to package
devtools::use_data(TDMaize, overwrite = TRUE)
