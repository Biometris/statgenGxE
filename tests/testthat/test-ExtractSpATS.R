context("Extract SpATS")

testTD <- createTD(data = testData[testData$field == "E1", ],
                   genotype = "seed", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoord = "Y", colCoord = "X")

modelSp <- STRunModel(testTD, design = "rowcol", traits = "t1")

test_that("the output of extract is of the proper type", {
  expect_is(STExtract(modelSp, what = "BLUEs"), "list")
  expect_is(STExtract(modelSp), "list")
  expect_length(STExtract(modelSp)[[1]], 17)
  expect_is(STExtract(modelSp, what = c("BLUEs", "BLUPs")), "list")
  expect_length(STExtract(modelSp, what = c("BLUEs", "BLUPs"))[[1]], 2)
})

extSp <- STExtract(modelSp)[[1]]
test_that("BLUEs are computed correctly", {
  expect_is(extSp$BLUEs, "data.frame")
  expect_identical(dim(extSp$BLUEs), c(15L, 2L))
  expect_equal(colnames(extSp$BLUEs), c("genotype", "t1"))
  expect_equal(extSp$BLUEs$t1,
               c(74.7151036023947, 84.6863159645867, 69.2058701507397,
                 72.4263349806975, 74.6775287748894, 84.3169610996821,
                 61.9016626250044, 76.4260688483351, 45.0089083719377,
                 55.73363594701, 127.386419996042, 35.3752613055704,
                 116.497846332665, 100.217877055618, 76.5100601280916))
})

test_that("SE of BLUEs are computed correctly", {
  expect_is(extSp$seBLUEs, "data.frame")
  expect_identical(dim(extSp$seBLUEs), c(15L, 2L))
  expect_equal(colnames(extSp$seBLUEs), c("genotype", "t1"))
  expect_equal(extSp$seBLUEs$t1,
               c(13.0937552217265, 12.7812399928054, 12.6750833951795,
                 12.9865301453826, 12.9605688922502, 13.0653342436885,
                 12.7352975349842, 12.6779389512144, 13.6714337388912,
                 13.126310040335, 15.3212096521888, 15.0833786817055,
                 12.4182933217468, 13.2613359160477, 12.3594075031988))
})

test_that("BLUPs are computed correctly", {
  expect_is(extSp$BLUPs, "data.frame")
  expect_identical(dim(extSp$BLUPs), c(15L, 2L))
  expect_equal(colnames(extSp$BLUPs), c("genotype", "t1"))
  expect_equal(extSp$BLUPs$t1,
               c(82.6087679270256, 76.8599245205723, 77.1948690330156,
                 82.4889022661131, 77.6036907961718, 77.4637164970361,
                 68.5420019338552, 75.5190085932238, 69.100508831864,
                 72.8694520315593, 93.4305094544051, 67.3168850301169,
                 97.2560475664543, 82.9993880603495, 76.1603123209831))
})

test_that("SE of BLUPs are computed correctly", {
  expect_is(extSp$seBLUPs, "data.frame")
  expect_identical(dim(extSp$seBLUPs), c(15L, 2L))
  expect_equal(colnames(extSp$seBLUPs), c("genotype", "t1"))
  expect_equal(extSp$seBLUPs$t1,
               c(10.9599851074304, 10.6712268311788, 10.9772900363026,
                 10.7832248077144, 11.0344062893833, 10.8689695115598,
                 10.8508399583286, 10.916541330236, 10.9347810564941,
                 10.7872768318315, 11.0977504304571, 10.8778327930184,
                 10.6709222386362, 10.6933088708978, 10.6547073855744))
})

test_that("heritability is computed correctly", {
  expect_is(extSp$heritability, "numeric")
  expect_length(extSp$heritability, 1)
  expect_equal(names(extSp$heritability), "t1")
  expect_equivalent(extSp$heritability, 0.46)
})

test_that("varGen is computed correctly", {
  expect_is(extSp$varGen, "numeric")
  expect_length(extSp$varGen, 1)
  expect_equal(names(extSp$varGen), "t1")
  expect_equivalent(extSp$varGen, 156.002441460485)
})

test_that("varSpat is computed correctly", {
  expect_is(extSp$varSpat, "matrix")
  expect_identical(dim(extSp$varSpat), c(5L, 1L))
  expect_equal(colnames(extSp$varSpat), "t1")
  expect_equal(as.numeric(extSp$varSpat),
               c(2022.23063560099, 136.171988031543,
                 7.50380433197137e-09, 0.000115160286180066,
                 0.117019622457697))
})

test_that("fitted values are computed correctly", {
  expect_is(extSp$fitted, "data.frame")
  expect_identical(dim(extSp$fitted), c(30L, 2L))
  expect_equal(colnames(extSp$fitted), c("genotype", "t1"))
  expect_equal(extSp$fitted$t1,
               c(71.2011919947091, 63.1039006172925, 55.7774636699599,
                 106.455101381074, 52.0527958005868, 49.7773596923711,
                 72.216949162496, 78.6978203767596, 91.4130323751844,
                 138.469400505856, 79.4720811382189, 58.3429716925874,
                 71.6927844745292, 77.2568018977556, 76.3260953095783,
                 85.6286934936668, 68.7505942521803, 60.7716588886156,
                 102.090734277441, 80.291466671846, 82.0370294498419,
                 70.6814124345326, 68.4586187576079, 81.0180757095049,
                 70.9985622813387, 49.8346571646361, 90.3661772058773,
                 105.805463114317, 111.185763325452, 59.4086442614029))
})

test_that("residuals are computed correctly", {
  expect_is(extSp$resid, "data.frame")
  expect_identical(dim(extSp$resid), c(30L, 2L))
  expect_equal(colnames(extSp$resid), c("genotype", "t1"))
  expect_equal(extSp$resid$t1,
               c(-6.20768167338727, 7.39043637551345, -16.224531377189,
                 -1.61750498202109, -13.7513442152902, 8.98003165678936,
                 13.7513442152901, 7.72427235215365, 3.78607568044147,
                 -3.78607568044146, -5.63937809977828, 5.63937809977823,
                 -7.39043637551346, 1.61750498202143, 16.6590847890083,
                 8.08123049403493, -6.51075834466288, -6.37487865373841,
                 16.224531377189, -8.98003165678931, -3.4260785879554,
                 3.42607858795543, -8.08123049403483, 6.51075834466296,
                 -16.6590847890082, 6.3748786537385, -7.72427235215368,
                 8.2108480989332, -8.2108480989332, 6.20768167338725))
})

test_that("rMeans are computed correctly", {
  expect_is(extSp$rMeans, "data.frame")
  expect_identical(dim(extSp$rMeans), c(30L, 2L))
  expect_equal(colnames(extSp$rMeans), c("genotype", "t1"))
  expect_equal(extSp$rMeans$t1,
               c(79.9923541999958, 60.29016783053, 56.7289964128857,
                 96.126879996754, 59.6937853822697, 56.6957894850249,
                 83.1663794313221, 75.9598539942815, 81.6957020300782,
                 111.056814747708, 73.1394524463974, 66.4380646346537,
                 85.6381282693448, 79.4423465994978, 69.3220286805539,
                 86.6221907998667, 74.6928561116224, 63.9984308682393,
                 92.2234803272238, 75.4125128764419, 74.3974673084901,
                 82.9399117320139, 73.3532106381913, 77.6473439897234,
                 81.2370204832617, 66.3036005968284, 85.1986832829286,
                 95.8588582258091, 91.5732633395968, 72.7377266556853))
})

test_that("random effects are computed correctly", {
  expect_is(extSp$ranEf, "data.frame")
  expect_identical(dim(extSp$ranEf), c(15L, 2L))
  expect_equal(colnames(extSp$ranEf), c("genotype", "t1"))
  expect_equal(extSp$ranEf$t1,
               c(4.11450226950924, -1.63434113694407, -1.29939662450081,
                 3.99463660859671, -0.890574861344598, -1.03054916048029,
                 -9.95226372366113, -2.9752570642926, -9.39375682565238,
                 -5.62481362595707, 14.9362437968887, -11.1773806273995,
                 18.7617819089379, 4.50512240283314, -2.33395333653327))
})

test_that("rDf is computed correctly", {
  expect_is(extSp$rDf, "numeric")
  expect_length(extSp$rDf, 1)
  expect_equal(names(extSp$rDf), "t1")
  expect_equivalent(extSp$rDf, 14)
})

test_that("effective dimensions are computed correctly", {
  expect_is(extSp$effDim, "matrix")
  expect_identical(dim(extSp$effDim), c(12L, 1L))
  expect_equal(colnames(extSp$effDim), "t1")
  expect_equal(as.numeric(extSp$effDim),
               c(1, 1, 1, 1, 6.40939170844998,
                 0.0688264628600227, 0.263826518026328,
                 0.27481928961094, 0.442528809466912,
                 9.93658502142113e-13, 1.9717150379426e-07,
                 1.31287578850783e-05))
})

test_that("ratios of effective dimensions are computed correctly", {
  expect_is(extSp$ratEffDim, "matrix")
  expect_identical(dim(extSp$ratEffDim), c(12L, 1L))
  expect_equal(colnames(extSp$ratEffDim), "t1")
  expect_equal(as.numeric(extSp$ratEffDim),
               c(1, 0.46, 0.01, 0.26, 1, 1, 1, 0.09, 0.07, 0, 0, 0))
})

test_that("correct attributes are added", {
  expect_equal(attr(x = extSp, which = "traits"), "t1")
  expect_equal(attr(x = extSp, which = "design"), "rowcol")
  expect_equal(attr(x = extSp, which = "engine"), "SpATS")
})
