context("Extract Lme4")

testTD <- createTD(data = testData[testData$field == "E1", ],
                   genotype = "seed", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoord = "Y", colCoord = "X")

modelLm <- STRunModel(testTD, design = "rcbd", traits = "t1", engine = "lme4")

test_that("the output of extract is of the proper type", {
  expect_is(STExtract(modelLm, what = "BLUEs"), "list")
  expect_is(STExtract(modelLm), "list")
  expect_length(STExtract(modelLm)[[1]], 18)
  expect_is(STExtract(modelLm, what = c("BLUEs", "BLUPs")), "list")
  expect_length(STExtract(modelLm, what = c("BLUEs", "BLUPs"))[[1]], 2)
})

extLm <- STExtract(modelLm)[[1]]
test_that("BLUEs are computed correctly", {
  expect_is(extLm$BLUEs, "data.frame"	)
  expect_identical(dim(extLm$BLUEs), c(15L, 2L))
  expect_equal(colnames(extLm$BLUEs), c("genotype", "t1"))
  expect_equal(extLm$BLUEs$t1,
               c(91.855951639415, 73.6623287954585, 74.8843349808426,
                 84.5319987913184, 68.9075264154031, 65.0344131821085,
                 55.3031580266259, 77.0436561256374, 62.1348724815412,
                 67.3983425459108, 108.495613219885, 65.304918128056,
                 114.94121644052, 78.9340989737003, 76.3592209421873))
})

test_that("SE of BLUEs are computed correctly", {
  expect_is(extLm$seBLUEs, "data.frame"	)
  expect_identical(dim(extLm$seBLUEs), c(15L, 2L))
  expect_equal(colnames(extLm$seBLUEs), c("genotype", "t1"))
  expect_equal(extLm$seBLUEs$t1, rep(x = 14.6535394925134, times = 15))
})

test_that("BLUPs are computed correctly", {
  expect_is(extLm$BLUPs, "data.frame"	)
  expect_identical(dim(extLm$BLUPs), c(15L, 2L))
  expect_equal(colnames(extLm$BLUPs), c("genotype", "t1"))
  expect_equal(extLm$BLUPs$t1,
               c(80.7705943535508, 76.7768113755654, 77.0450607255733,
                 79.1628728751421, 75.7330583544096, 74.8828498041054,
                 72.746688055138, 77.5190652768411, 74.2463555559162,
                 75.4017690394489, 84.4232586826672, 74.9422298439057,
                 85.838168766613, 77.9340468630927, 77.3688211166405))
})

test_that("SE of BLUPs are computed correctly", {
  expect_is(extLm$seBLUPs, "data.frame"	)
  expect_identical(dim(extLm$seBLUPs), c(15L, 2L))
  expect_equal(colnames(extLm$seBLUPs), c("genotype", "t1"))
  expect_equal(extLm$seBLUPs$t1, rep(x = 6.86554949586288, times = 15))
})

test_that("unit errors are computed correctly", {
  expect_is(extLm$ue, "data.frame"	)
  expect_identical(dim(extLm$ue), c(15L, 2L))
  expect_equal(colnames(extLm$ue), c("genotype", "t1"))
  expect_equal(extLm$ue$t1, rep(x = 214.726219658649, times = 15))
})

test_that("heritability is computed correctly", {
  expect_is(extLm$heritability, "numeric")
  expect_length(extLm$heritability, 1)
  expect_equal(names(extLm$heritability), "t1")
  expect_equivalent(extLm$heritability, 0.219515541914842)
})

test_that("varGen is computed correctly", {
  expect_is(extLm$varGen, "numeric")
  expect_length(extLm$varGen, 1)
  expect_equal(names(extLm$varGen), "t1")
  expect_equivalent(extLm$varGen, 60.3929641286978)
})

test_that("varErr is computed correctly", {
  expect_is(extLm$varErr, "numeric")
  expect_length(extLm$varErr, 1)
  expect_equal(names(extLm$varErr), "t1")
  expect_equivalent(extLm$varErr, 429.452689034918)
})

test_that("fitted values are computed correctly", {
  expect_is(extLm$fitted, "data.frame"	)
  expect_identical(dim(extLm$fitted), c(30L, 3L))
  expect_equal(colnames(extLm$fitted), c("genotype", "repId", "t1"))
  expect_equal(extLm$fitted$t1,
               c(98.1557845869273, 83.3434890731497, 68.4347054290536,
                 73.6981754934231, 114.795446167397, 71.6047510755682,
                 121.241049388032, 85.2339319212126, 82.6590538896995,
                 79.9621617429708, 81.1841679283549, 90.8318317388307,
                 75.2073593629154, 71.3342461296208, 61.6029909741381,
                 85.5561186919027, 70.7438231781252, 55.835039534029,
                 61.0985095983986, 102.195780272372, 59.0050851805438,
                 108.641383493008, 72.6342660261881, 70.059387994675,
                 67.3624958479463, 68.5845020333304, 78.2321658438062,
                 62.6076934678909, 58.7345802345963, 49.0033250791136))
})

test_that("residuals are computed correctly", {
  expect_is(extLm$resid, "data.frame"	)
  expect_identical(dim(extLm$resid), c(30L, 3L))
  expect_equal(colnames(extLm$resid), c("genotype", "repId", "t1"))
  expect_equal(extLm$resid$t1,
               c(6.68181181212573, 10.3664349145521, 17.5335879487325,
                 -9.39582739440736, -0.779134954146725, -6.61124075424643,
                 13.442275437382, 33.0813337334172, -8.55156286721146,
                 -25.6226842506403, -18.9443320208374, -4.40973900991743,
                 -1.37465632447471, -0.0228111145641388, -5.39345515576354,
                 -6.68181181212574, -10.3664349145521, -17.5335879487324,
                 9.39582739440736, 0.779134954146727, 6.61124075424643,
                 -13.442275437382, -33.0813337334172, 8.55156286721146,
                 25.6226842506403, 18.9443320208374, 4.40973900991743,
                 1.37465632447471, 0.0228111145641317, 5.39345515576353))
})

test_that("standardized residuals are computed correctly", {
  expect_is(extLm$stdRes, "data.frame"	)
  expect_identical(dim(extLm$stdRes), c(30L, 3L))
  expect_equal(colnames(extLm$stdRes), c("genotype", "repId", "t1"))
  expect_equal(extLm$stdRes$t1,
               c(0.4719905598143, 0.732265372951495, 1.23853951954428,
                 -0.663703491881794, -0.0550366208325776, -0.467005553643564,
                 0.949536935080114, 2.33680290125858, -0.604066241084128,
                 -1.80993799637666, -1.33819181492009, -0.311495630593602,
                 -0.0971031250780877, -0.00161133402673469, -0.380983480211527,
                 -0.4719905598143, -0.732265372951495, -1.23853951954428,
                 0.663703491881794, 0.0550366208325777, 0.467005553643564,
                 -0.949536935080114, -2.33680290125858, 0.604066241084128,
                 1.80993799637666, 1.33819181492009, 0.311495630593602,
                 0.0971031250780872, 0.00161133402673419, 0.380983480211527))
})

test_that("rMeans are computed correctly", {
  expect_is(extLm$rMeans, "data.frame"	)
  expect_identical(dim(extLm$rMeans), c(30L, 3L))
  expect_equal(colnames(extLm$rMeans), c("genotype", "repId", "t1"))
  expect_equal(extLm$rMeans$t1,
               c(87.0704273010631, 83.8188982243534, 80.5461885034285,
                 81.7016019869611, 90.7230916301794, 81.242062791418,
                 92.1380017141253, 84.233879810605, 83.6686540641528,
                 83.0766443230777, 83.3448936730856, 85.4627058226544,
                 82.0328913019219, 81.1826827516177, 79.0465210026503,
                 74.4707614060385, 71.2192323293289, 67.946522608404,
                 69.1019360919366, 78.1234257351549, 68.6423968963934,
                 79.5383358191008, 71.6342139155804, 71.0689881691282,
                 70.4769784280532, 70.7452277780611, 72.8630399276299,
                 69.4332254068973, 68.5830168565931, 66.4468551076258))
})

test_that("random effects are computed correctly", {
  expect_is(extLm$ranEf, "data.frame"	)
  expect_identical(dim(extLm$ranEf), c(15L, 2L))
  expect_equal(colnames(extLm$ranEf), c("genotype", "t1"))
  expect_equal(extLm$ranEf$t1,
               c(3.1178176409768, -0.875965337008561, -0.60771598700067,
                 1.51009616256813, -1.9197183581644, -2.76992690846861,
                 -4.90608865743599, -0.13371143573287, -3.40642115665777,
                 -2.25100767312514, 6.77048197009317, -2.71054686866831,
                 8.18539205403902, 0.281270150518691, -0.283955595933513))
})

test_that("Waldtest is computed correctly", {
  expect_is(extLm$wald$t1, "data.frame")
  expect_length(extLm$wald$t1, 4)
  expect_equivalent(unlist(extLm$wald$t1),
                    c(15, 14, 29.278, 5.50452147204485e-08))
})

test_that("CV is computed correctly", {
  expect_is(extLm$CV, "numeric")
  expect_length(extLm$CV, 1)
  expect_equal(names(extLm$CV), "t1")
  expect_equivalent(extLm$CV, 26.6870486341882)
})

test_that("rDf is computed correctly", {
  expect_is(extLm$rDf, "integer")
  expect_length(extLm$rDf, 1)
  expect_equal(names(extLm$rDf), "t1")
  expect_equivalent(extLm$rDf, 14)
})

test_that("correct attributes are added", {
  expect_equal(attr(x = extLm, which = "traits"), "t1")
  expect_equal(attr(x = extLm, which = "design"), "rcbd")
  expect_equal(attr(x = extLm, which = "engine"), "lme4")
})

