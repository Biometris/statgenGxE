context("Modeling")

testTD <- createTD(data = testData[testData$field == "E1", ],
                   genotype = "seed", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoordinates = "Y", colCoordinates = "X")

modelSp <- STRunModel(testTD, design = "rowcol", traits = "t1")
test_that("running models creates objects of proper class - SpATS", {
  expect_is(modelSp, "SSA")
  expect_is(modelSp$mFix$t1, "SpATS")
  expect_is(modelSp$mRand$t1, "SpATS")
  expect_equal(modelSp$traits, "t1")
  expect_equal(modelSp$design, "rowcol")
  expect_equal(modelSp$engine, "SpATS")
})

modelLm <- STRunModel(testTD, design = "rcbd", traits = "t1")
test_that("running models creates objects of proper class - lme4", {
  expect_is(modelLm$mFix$t1, "lm")
  expect_is(modelLm$mRand$t1, "lmerMod")
  expect_equal(modelLm$design, "rcbd")
  expect_equal(modelLm$engine, "lme4")
})

modelAs <- STRunModel(testTD, design = "res.ibd", traits = "t1", engine = "asreml")
test_that("running models creates objects of proper class - asreml", {
  skip_on_cran()
  expect_is(modelAs$mFix$t1, "asreml")
  expect_is(modelAs$mRand$t1, "asreml")
  expect_equal(modelAs$design, "res.ibd")
  expect_equal(modelAs$engine, "asreml")
})

modelSpF <- STRunModel(testTD, design = "rowcol", traits = "t1", what = "fixed")
modelSpR <- STRunModel(testTD, design = "rowcol", traits = "t1", what = "random")
test_that("option what produces expected output", {
  expect_null(modelSpF$mRand)
  expect_is(modelSpF$mFix$t1, "SpATS")
  expect_null(modelSpR$mFix)
  expect_is(modelSpR$mRand$t1, "SpATS")
})


