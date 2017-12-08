## Test SES 2016 for field ARC601
rm(list = ls())
## Load cleaned data for analysis.
load("D:/R packages/SESvdH/data/2016/cleandata.RData")
rm(list = setdiff(ls(), c("dfARC601", "traitnames")))
## Load results for comparison
load("D:/R packages/SESvdH/data/2016/BLUEs_per_field.RData")
load("D:/R packages/SESvdH/data/2016/RCcorrected_per_field.RData")
load("D:/R packages/SESvdH/results/2016/heritabilities.RData")

## Create TD object.
TDARC601 <- createTD(data = dfARC601, genotype = "Seedname",
                     rowId = "Yf", colId = "Xf", rowCoordinates = "Y",
                     colCoordinates = "X", design = "rowcol")

## Perform SpATS analysis for computing BLUEs. - step 5 in SES workflow
ARC601Sp <- STRunModel(TDARC601, traits = traitnames, engine = "SpATS")
BLUEsSp <- STExtract(ARC601Sp, what = "BLUEs")

## Comparison for some traits.
## Small differences because of difference in number of segments
## for rows in SpATS analysis.
BLUEsOrig <- BLUEs[BLUEs$FieldId == "ARC601", ]
plot(x = BLUEsSp$THA, y = BLUEsOrig$THA)
abline(a = 0, b = 1, col = "blue")
plot(x = BLUEsSp$pctS, y = BLUEsOrig$pctS)
abline(a = 0, b = 1, col = "blue")

## Perform analysis using asreml and lme4. - step 7 in SES workflow
## lme4 used for extra comparison.
## traits restricted to "THA" for speed
ARC601As <- STRunModel(TDARC601, traits = "THA", engine = "asreml")
BLUEsAs <- STExtract(ARC601As, what = "BLUEs")
ARC601Lm <- STRunModel(TDARC601, traits = "THA", engine = "lme4")
BLUEsLm <- STExtract(ARC601Lm, what = "BLUEs")

## Comparison for some THA.
## Plots seem perfect, but there are small differences when using asreml.
## This is because in SES the model is fitted in one step whereas in RAP
## it is done in two steps and some variance components are fixed in between.
RCCorrectedOrig <- RCcorrected[RCcorrected$FieldId == "ARC601", ]
meanDifAs <- mean(abs(RCCorrectedOrig$THA - BLUEsAs$THA))
meanDifLm <- mean(abs(RCCorrectedOrig$THA - BLUEsLm$THA))
plot(x = BLUEsAs$THA, y = RCCorrectedOrig$THA)
abline(a = 0, b = 1, col = "blue")
plot(x = BLUEsLm$THA, y = RCCorrectedOrig$THA)
abline(a = 0, b = 1, col = "blue")

## Comparison for heritabilities. - step 10 in SES workflow
## Again small differences because of difference in number of segments
## for rows in SpATS analysis.
herit <- STExtract(ARC601Sp, what = "heritability")
heritOrig <- heritabilities["ARC601", ]
plot(x = herit, y = heritOrig)
abline(a = 0, b = 1, col = "blue")

