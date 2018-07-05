context("Class SSA")

## Create testdata containing only one trial.
testTD <- createTD(data = testData, trial = "field",
                   genotype = "seed", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoord = "Y", colCoord = "X")
modelSp <- STRunModel(testTD, trials = "E1", design = "rowcol", traits = "t1")
test_that("summary.SSA produces correct output for SpATS", {
  sumSp <- summary(modelSp)
  expect_length(sumSp, 6)
  expect_null(sumSp$selSpatMod)
  expect_equal(nrow(sumSp$stats), 9)
  expect_equal(dim(sumSp$meanTab), c(15, 4))
  expect_equal(unname(sumSp$heritability), 0.46)
  expect_equal(nrow(sumSp$sed), 0)
  expect_equal(nrow(sumSp$lsd), 0)
})

test_that("summary.SSA produces correct output for lme4", {
  modelLm <- STRunModel(testTD, trials = "E1", design = "rowcol", traits = "t1",
                        engine = "lme4")
  sumLm <- summary(modelLm)
  expect_length(sumLm, 6)
  expect_null(sumLm$selSpatMod)
  expect_equal(nrow(sumLm$stats), 9)
  expect_equal(dim(sumLm$meanTab), c(15, 4))
  expect_equal(unname(sumLm$heritability), 0.121957305492988)
  expect_equal(nrow(sumLm$sed), 0)
  expect_equal(nrow(sumLm$lsd), 0)
})

test_that("summary.SSA produces correct output for asreml", {
  skip_on_cran()
  modelAs <- STRunModel(testTD, trials = "E1", design = "rowcol", traits = "t1",
                        engine = "asreml")
  sumAs <- summary(modelAs)
  expect_length(sumAs, 6)
  expect_null(sumAs$selSpatMod)
  expect_equal(nrow(sumAs$stats), 9)
  expect_equal(dim(sumAs$meanTab), c(15, 4))
  expect_equal(unname(sumAs$heritability), 0.206091659337936)
  expect_equal(nrow(sumAs$sed), 3)
  expect_equal(nrow(sumAs$lsd), 3)
})

test_that("option nBest functions properly", {
  sumSp <- summary(modelSp, nBest = 5)
  expect_equal(dim(sumSp$meanTab), c(5, 4))
})

