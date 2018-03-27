context("Modeling")

## Create testdata containing only one trial.
testTD <- createTD(data = testData, trial = "field",
                   genotype = "seed", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoordinates = "Y", colCoordinates = "X")

## Helper function for testing base structure that has to be consistent
## for all SSA objects independent of engine and options.
expect_SSA <- function(SSA) {
  test_that(paste(deparse(substitute(SSA)), "has correct SSA structure"), {
    expect_is(SSA, "SSA")
    for (tr in names(SSA)) {
      expect_length(SSA[[tr]], 8)
      expect_named(SSA[[tr]], c("mRand", "mFix", "TD", "traits", "design",
                                "spatial", "engine", "predicted"))
      expect_is(SSA[[tr]]$TD, "TD")
    }
  })
}

## Helper function for testing base structure for fitted models within an
## SSA object. Param class is used to indicate the output type of the fitted
## model. Normally this is identical to the engine in SSA but for lme4 this
## can vary depending on the fitted model.
expect_SSAMod <- function(SSA,
                          what,
                          class = NULL) {
  if (is.null(class)) {
    class <- SSA[[1]]$engine
  }
  for (tr in names(SSA)) {
    SSAMod <- SSA[[tr]][[what]]
    test_that(paste(deparse(substitute(what)), "in", deparse(substitute(SSA)),
                    "has correct structure"), {
                      expect_is(SSAMod, "list")
                      expect_length(SSAMod, length(SSA[[tr]]$traits))
                      expect_named(SSAMod, SSA[[tr]]$traits)
                      for (trait in SSA[[tr]]$traits) {
                        expect_is(SSAMod[[trait]], class)
                      }
                    })
  }
}

test_that("running models creates objects with correct structure - SpATS", {
  modelSp <- STRunModel(testTD, trials = "E1", design = "rowcol", traits = "t1")
  expect_SSA(modelSp)
  expect_SSAMod(modelSp, "mRand")
  expect_SSAMod(modelSp, "mFix")
  expect_identical(modelSp[["E1"]]$traits, "t1")
  expect_identical(modelSp[["E1"]]$design, "rowcol")
  expect_identical(modelSp[["E1"]]$spatial, "2 dimensional P-splines")
  expect_identical(modelSp[["E1"]]$engine, "SpATS")
})

test_that("running models creates objects with correct structure - lme4", {
  modelLm <- STRunModel(testTD, trials = "E1", design = "rcbd", traits = "t1")
  expect_SSA(modelLm)
  expect_SSAMod(modelLm, "mRand", class = "lmerMod")
  expect_SSAMod(modelLm, "mFix", class = "lm")
  expect_identical(modelLm[["E1"]]$traits, "t1")
  expect_identical(modelLm[["E1"]]$design, "rcbd")
  expect_false(modelLm[["E1"]]$spatial)
  expect_identical(modelLm[["E1"]]$engine, "lme4")
})

test_that("running models creates objects with correct structure - asreml", {
  skip_on_cran()
  modelAs <- STRunModel(testTD, trials = "E1", design = "res.ibd",
                        traits = "t1", engine = "asreml")
  expect_SSA(modelAs)
  expect_SSAMod(modelAs, "mRand")
  expect_SSAMod(modelAs, "mFix")
  expect_identical(modelAs[["E1"]]$traits, "t1")
  expect_identical(modelAs[["E1"]]$design, "res.ibd")
  expect_false(modelAs[["E1"]]$spatial)
  expect_identical(modelAs[["E1"]]$engine, "asreml")
})

test_that("option what produces expected output - SpATS", {
  modelSp <- STRunModel(testTD, trials = "E1", design = "rowcol", traits = "t1")
  modelSpF <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = "t1", what = "fixed")
  expect_SSA(modelSpF)
  expect_null(modelSpF[["E1"]]$mRand)
  expect_SSAMod(modelSpF, "mFix")
  expect_equal(modelSpF[["E1"]]$mFix, modelSp[["E1"]]$mFix)
  modelSpR <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = "t1", what = "random")
  expect_SSA(modelSpR)
  expect_SSAMod(modelSpR, "mRand")
  expect_equal(modelSp[["E1"]]$mRand, modelSpR[["E1"]]$mRand)
  expect_null(modelSpR[["E1"]]$mFix)
})

test_that("option what produces expected output - lme4", {
  modelLm <- STRunModel(testTD, trials = "E1", design = "rcbd", traits = "t1")
  modelLmF <- STRunModel(testTD, trials = "E1", design = "rcbd",
                         traits = "t1", what = "fixed")
  expect_SSA(modelLmF)
  expect_null(modelLmF[["E1"]]$mRand)
  expect_SSAMod(modelLmF, "mFix", "lm")
  expect_equal(modelLmF[["E1"]]$mFix, modelLm[["E1"]]$mFix)
  modelLmR <- STRunModel(testTD, trials = "E1", design = "rcbd",
                         traits = "t1", what = "random")
  expect_SSA(modelLmR)
  expect_SSAMod(modelLmR, "mRand", "lmerMod")
  expect_equal(modelLm[["E1"]]$mRand, modelLmR[["E1"]]$mRand)
  expect_null(modelLmR[["E1"]]$mFix)
})

test_that("option what produces expected output - asreml", {
  skip_on_cran()
  modelAs <- STRunModel(testTD, trials = "E1", design = "rowcol", traits = "t1",
                        engine = "asreml")
  modelAsF <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = "t1", what = "fixed", engine = "asreml")
  expect_SSA(modelAsF)
  expect_null(modelAsF[["E1"]]$mRand)
  expect_SSAMod(modelAsF, "mFix")
  modelAsR <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = "t1", what = "random", engine = "asreml")
  expect_SSA(modelAsR)
  expect_SSAMod(modelAsR, "mRand")
  expect_equal(modelAs[["E1"]]$mRand, modelAsR[["E1"]]$mRand)
  expect_null(modelAsR[["E1"]]$mFix)
})

test_that("running models for multiple traits produces correct output structure", {
  modelSp2 <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = paste0("t", 1:2))
  modelSp3 <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = paste0("t", 1:3))
  modelSp4 <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = paste0("t", 1:4))
  expect_SSA(modelSp2)
  expect_SSA(modelSp3)
  expect_SSA(modelSp4)
  expect_SSAMod(modelSp2, "mRand")
  expect_SSAMod(modelSp3, "mRand")
  expect_SSAMod(modelSp4, "mRand")
  expect_SSAMod(modelSp2, "mFix")
  expect_SSAMod(modelSp3, "mFix")
  expect_SSAMod(modelSp4, "mFix")
})

test_that("running models for multiple traits doesn't change trait results", {
  modelSp1 <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = "t1")
  modelSp2 <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = paste0("t", 1:2))
  modelSp3 <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = paste0("t", 1:3))
  modelSp4 <- STRunModel(testTD, trials = "E1", design = "rowcol",
                         traits = paste0("t", 1:4))
  expect_equal(modelSp1[["E1"]]$mRand, modelSp2[["E1"]]$mRand["t1"])
  expect_equal(modelSp2[["E1"]]$mRand, modelSp3[["E1"]]$mRand[paste0("t", 1:2)])
  expect_equal(modelSp3[["E1"]]$mRand, modelSp4[["E1"]]$mRand[paste0("t", 1:3)])
  expect_equal(modelSp1[["E1"]]$mFix, modelSp2[["E1"]]$mFix["t1"])
  expect_equal(modelSp2[["E1"]]$mFix, modelSp3[["E1"]]$mFix[paste0("t", 1:2)])
  expect_equal(modelSp3[["E1"]]$mFix, modelSp4[["E1"]]$mFix[paste0("t", 1:3)])
})

test_that("option covariates produces expected output structure", {
  modelSpCov <- STRunModel(testTD, trials = "E1", design = "rowcol",
                           traits = "t1", covariates = "repId")
  expect_SSA(modelSpCov)
  expect_SSAMod(modelSpCov, "mRand")
  expect_SSAMod(modelSpCov, "mFix")
  expect_gt(grep(pattern = "repId",
                 x = deparse(modelSpCov[["E1"]]$mRand$t1$model$fixed)), 0)
  expect_gt(grep(pattern = "repId",
                 x = deparse(modelSpCov[["E1"]]$mFix$t1$model$fixed)), 0)
  modelLmCov <- STRunModel(testTD, trials = "E1", design = "rcbd",
                           traits = "t1", covariates = "repId")
  expect_SSA(modelLmCov)
  expect_SSAMod(modelLmCov, "mRand", "lmerMod")
  expect_SSAMod(modelLmCov, "mFix", "lm")
  expect_true("repId" %in% colnames(modelLmCov[["E1"]]$mRand$t1@frame))
  expect_true("repId" %in% colnames(modelLmCov[["E1"]]$mFix$t1$model))
  skip_on_cran()
  modelAsCov <- STRunModel(testTD, trials = "E1", design = "rowcol",
                           traits = "t1", covariates = "repId",
                           engine = "asreml")
  expect_SSA(modelAsCov)
  expect_SSAMod(modelAsCov, "mRand")
  expect_SSAMod(modelAsCov, "mFix")
  expect_gt(grep(pattern = "repId",
                 x = deparse(modelAsCov[["E1"]]$mRand$t1$fixed.formula)), 0)
  expect_gt(grep(pattern = "repId",
                 x = deparse(modelAsCov[["E1"]]$mFix$t1$fixed.formula)), 0)
})

test_that("option useCheckId produces expected output structure", {
  modelSpCi <- STRunModel(testTD, trials = "E1", design = "rowcol", traits = "t1",
                          useCheckId = TRUE)
  expect_SSA(modelSpCi)
  expect_SSAMod(modelSpCi, "mRand")
  expect_SSAMod(modelSpCi, "mFix")
  expect_gt(grep(pattern = "checkId",
                 x = deparse(modelSpCi[["E1"]]$mRand$t1$model$fixed)), 0)
  expect_gt(grep(pattern = "checkId",
                 x = deparse(modelSpCi[["E1"]]$mFix$t1$model$fixed)), 0)
  modelLmCi <- STRunModel(testTD, trials = "E1", design = "rcbd", traits = "t1",
                          useCheckId = TRUE)
  expect_SSA(modelLmCi)
  expect_SSAMod(modelLmCi, "mRand", "lmerMod")
  expect_SSAMod(modelLmCi, "mFix", "lm")
  expect_true("checkId" %in% colnames(modelLmCi[["E1"]]$mRand$t1@frame))
  expect_true("checkId" %in% colnames(modelLmCi[["E1"]]$mFix$t1$model))
})

test_that("option trySpatial produces expected output structure", {
  modelSp <- STRunModel(testTD, trials = "E1", design = "rowcol", traits = "t1")
  modelSpTs <- STRunModel(testTD, trials = "E1", design = "rowcol",
                          traits = "t1", trySpatial = TRUE)
  expect_SSA(modelSpTs)
  expect_SSAMod(modelSpTs, "mRand")
  expect_SSAMod(modelSpTs, "mFix")
  ## SpATS should use trySpatial as default. Timestamp will be different.
  expect_equivalent(modelSp, modelSpTs)
  expect_warning(STRunModel(testTD, trials = "E1", design = "rowcol",
                            traits = "t1", trySpatial = TRUE, engine = "lme4"),
                 "Spatial models can only be fitted using SpATS or asreml.")
  skip_on_cran()
  modelAsTs <- STRunModel(testTD, trials = "E1", design = "rowcol",
                          traits = "t1", trySpatial = TRUE, engine = "asreml")
  expect_SSA(modelAsTs)
  expect_SSAMod(modelAsTs, "mRand")
  expect_SSAMod(modelAsTs, "mFix")
  expect_is(modelAsTs[["E1"]]$spatial, "list")
  expect_length(modelAsTs[["E1"]]$spatial, 1)
  expect_named(modelAsTs[["E1"]]$spatial, "t1")
  expect_identical(modelAsTs[["E1"]]$spatial$t1, "id(x)exp")
})

