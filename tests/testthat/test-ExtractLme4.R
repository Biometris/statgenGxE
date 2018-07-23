context("Extract Lme4")

testTD <- createTD(data = testData[testData$field == "E1", ],
                   genotype = "seed", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoord = "Y", colCoord = "X")

modelLm <- STRunModel(testTD, design = "rcbd", traits = "t1")

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
               c(80.770606359766, 76.7768080023636, 77.0450583853567,
                 79.16287869028, 75.7330509618824, 74.8828391375616,
                 72.7466691625783, 77.5190647619399, 74.2463424383353,
                 75.4017603711797, 84.4232847547063, 74.9422194060249,
                 85.8382002872437, 77.9340479462189, 77.3688200231731))
})

test_that("SE of BLUPs are computed correctly", {
  expect_is(extLm$seBLUPs, "data.frame"	)
  expect_identical(dim(extLm$seBLUPs), c(15L, 2L))
  expect_equal(colnames(extLm$seBLUPs), c("genotype", "t1"))
  expect_equal(extLm$seBLUPs$t1, rep(x = 6.86556085594243, times = 15))
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
  expect_equivalent(extLm$heritability, 0.219516387233952)
})

test_that("varGen is computed correctly", {
  expect_is(extLm$varGen, "numeric")
  expect_length(extLm$varGen, 1)
  expect_equal(names(extLm$varGen), "t1")
  expect_equivalent(extLm$varGen, 60.3932293973456)
})

test_that("varErr is computed correctly", {
  expect_is(extLm$varErr, "numeric")
  expect_length(extLm$varErr, 1)
  expect_equal(names(extLm$varErr), "t1")
  expect_equivalent(extLm$varErr, 429.452456471173)
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
               c(87.0704393072783, 83.8188977094521, 80.5461753858475,
                 81.701593318692,  90.7231177022186, 81.2420523535372,
                 92.138033234756, 84.2338808937312, 83.6686529706854,
                 83.0766409498759, 83.3448913328689, 85.4627116377923,
                 82.0328839093947, 81.1826720850739, 79.0465021100905,
                 74.4707734122537, 71.2192318144276, 67.946509490823,
                 69.1019274236675, 78.123451807194, 68.6423864585126,
                 79.5383673397315, 71.6342149987067, 71.0689870756608,
                 70.4769750548513, 70.7452254378444, 72.8630457427678,
                 69.4332180143702, 68.5830061900493, 66.446836215066))
})

test_that("random effects are computed correctly", {
  expect_is(extLm$ranEf, "data.frame"	)
  expect_identical(dim(extLm$ranEf), c(15L, 2L))
  expect_equal(colnames(extLm$ranEf), c("genotype", "t1"))
  expect_equal(extLm$ranEf$t1,
               c(3.11782964719197, -0.875968710210448, -0.607718327217375,
                 1.510101977706, -1.9197257506916, -2.76993757501244,
                 -4.90610754999576, -0.133711950634149, -3.40643427423876,
                 -2.25101634139431, 6.77050804213227, -2.71055730654913,
                 8.18542357466969, 0.28127123364488, -0.283956689400929))
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

