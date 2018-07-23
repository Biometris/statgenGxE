library(RAP)
## Get latest version
devtools::install_github("CIP-RIU/brapi")

## BrAPI Connection object for Authentication on BMS
bmscon <- brapi::ba_connect(brapiDb = NULL, secure = FALSE, protocol = "http://",
                            db = "34.226.132.187", port = 48080,
                            apipath = "bmsapi", multicrop = TRUE, crop = "",
                            user = "wur", password = "27Aqua74", token = "",
                            granttype = "password", clientid = "rbrapi",
                            bms = TRUE)
## Authenticate user with password and obtain the session token
bmscon <- brapi::ba_login(con = bmscon)
## Select "wheat" as crop and write to connection object
bmscon$crop <- "wheat"
## Extract data for all fields (studyDBIDs "21" to "38")
studyDbIds <- as.character(21:38)  #  38)
studyTablesList <- lapply(studyDbIds, function(id) {
  brapi::ba_studies_table(con = bmscon, studyDbId = id)
})
## Create one data.frame with data for all fields.
studyTableTot <- Reduce("rbind", studyTablesList)
## Remove data for sites 3, 13, 15 and 17 (failed crop)
studyTableTot <- studyTableTot[!studyTableTot$locationName %in%
                                 c("Site03", "Site13", "Site15", "Site17"), ]
studyTableTot$locationName <- droplevels(studyTableTot$locationName)
## Rename trait column for ease of use
colnames(studyTableTot)[colnames(studyTableTot) == "GY_Calc_tha|22661"] <-
  "GY_Calc_tha"
## Create TD object
wheatTD <- createTD(data = studyTableTot, genotype = "germplasmName",
                    trial = "locationName", repId = "replicate",
                    subBlock = "blockNumber", trDesign = "rcbd")
### Single Site Analysis
## Run Single Site Analysis using asreml
SSA <- STRunModel(TD = wheatTD, traits = "GY_Calc_tha", engine = "asreml")
## Create a report for the single site analysis.
report(SSA, outfile = "./testReports/SSASite01.pdf", trial = "Site01")
## Create new TD objects from BLUPs and BLUEs
wheatTDBlup <- SSAtoTD(SSA, what = "BLUPs")
wheatTDBlue <- SSAtoTD(SSA, what = c("BLUEs", "seBLUEs"), addWt = TRUE)

###GxE analysis
## Run AMMI analysis
ammi <- gxeAmmi(TD = wheatTDBlup, trait = "GY_Calc_tha")
summary(ammi)
report(ammi, outfile = "./testReports/ammiWheat.pdf")

## Run Finlay-Wilkinson analysis
fw <- gxeFw(TD = wheatTDBlup, trait = "GY_Calc_tha", sorted = "ascending")
summary(fw)
report(fw, outfile = "./testReports/fwWheat.pdf")

## Run stability analysis
stab <- gxeStability(TD = wheatTDBlup, trait = "GY_Calc_tha")
summary(stab)
report(stab, outfile = "./testReports/stabWheat.pdf")

## Run variance covariance analysis to find the best covariance structure.
varComp <- gxeVarComp(TD = wheatTDBlup, trait = "GY_Calc_tha",
                      engine = "asreml")
summary(varComp)
varComp2 <- gxeVarComp(TD = wheatTDBlue, trait = "BLUEs_GY_Calc_tha",
                       engine = "asreml", useWt = TRUE)
summary(varComp2)
report(varComp, outfile = "./testReports/varCompWheat.pdf")
