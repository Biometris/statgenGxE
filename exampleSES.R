## Test SES 2016 for field ARC601 and BAU601
rm(list = ls())
## Load cleaned data for analysis.
load("D:/R packages/SESvdH/data/2016/cleandata.RData")
rm(list = setdiff(ls(), c("dfARC601", "dfBAU601", "traitnames")))
## Load results for comparison
load("D:/R packages/SESvdH/data/2016/BLUEs_per_field.RData")
load("D:/R packages/SESvdH/data/2016/RCcorrected_per_field.RData")
load("D:/R packages/SESvdH/results/2016/heritabilities.RData")

## Create TD object.
TDARC601 <- createTD(data = rbind(dfARC601, dfBAU601),
                     genotype = "Seedname", trial = "FieldId",
                     rowCoord = "Y", colCoord = "X", trDesign = "rowcol")
plot(TDARC601)

## Perform SpATS analysis for computing BLUEs. - step 5 in SES workflow
ARC601Sp <- STRunModel(TDARC601, traits = "THA", engine = "SpATS",
                       control = list(nSeg = c(ceiling(nlevels(TDARC601$ARC601$colId) / 2),
                                               nlevels(TDARC601$ARC601$rowId))))
BLUEsSp <- STExtract(ARC601Sp, what = "BLUEs")$ARC601$BLUEs

## Comparison for THA and pctS.
BLUEsOrig <- BLUEs[BLUEs$FieldId == "ARC601", ]
plot(x = BLUEsSp$THA, y = BLUEsOrig$THA)
abline(a = 0, b = 1, col = "blue")
#plot(x = BLUEsSp$pctS, y = BLUEsOrig$pctS)
#abline(a = 0, b = 1, col = "blue")

## Perform analysis using asreml and lme4. - step 7 in SES workflow
## traits restricted to "THA" for speed
## For asreml model only fitted with genotype fixed replicate results
ARC601As <- STRunModel(TDARC601, traits = "THA", engine = "asreml", what = "fixed")
BLUEsAs <- STExtract(ARC601As, what = "BLUEs")$ARC601$BLUEs
## Model also fitted using for extra comparison.
ARC601Lm <- STRunModel(TDARC601, traits = "THA", engine = "lme4")
BLUEsLm <- STExtract(ARC601Lm, what = "BLUEs")$ARC601$BLUEs

## Comparison for THA.
RCCorrectedOrig <- RCcorrected[RCcorrected$FieldId == "ARC601", ]
plot(x = BLUEsAs$THA, y = RCCorrectedOrig$THA)
abline(a = 0, b = 1, col = "blue")
plot(x = BLUEsLm$THA, y = RCCorrectedOrig$THA)
abline(a = 0, b = 1, col = "blue")

## Comparison for heritabilities. - step 10 in SES workflow
## Create new TD object since for this calculation column Geno is used as genotype
## and "Check" as a fixed factor
TDARC601H <- createTD(data = dfARC601, genotype = "Geno", checkId = "Check",
                     rowCoord = "Y", colCoord = "X", trDesign = "rowcol")
## Only fit model with genotype random.
## Because of the structure of the data, genotype fixed leads to collinearity
ARC601SpH <- STRunModel(TDARC601H, traits = traitnames, what = "random",
                        useCheckId = TRUE, engine = "SpATS",
                        control = list(nSeg = c(ceiling(nlevels(TDARC601H$dfARC601$colId) / 2),
                                                nlevels(TDARC601H$dfARC601$rowId))))
herit <- STExtract(ARC601SpH, what = "heritability")$dfARC601$heritability
heritOrig <- heritabilities["ARC601", ]
plot(x = herit, y = heritOrig)
abline(a = 0, b = 1, col = "blue")


