context("Modeling")

## Create testdata containing only one trial.
testTD <- createTD(data = testData, trial = "field", genotype = "seed",
                   repId = "rep", subBlock = "block", rowCoord = "Y",
                   colCoord = "X")

## Helper function for testing base structure that has to be consistent
## for all SSA objects independent of engine and options.
expect_SSA <- function(SSA) {
  test_that(paste(deparse(substitute(SSA)), "has correct SSA structure"), {
    expect_is(SSA, "SSA")
    for (tr in names(SSA)) {
      expect_length(SSA[[tr]], 9)
      expect_named(SSA[[tr]], c("mRand", "mFix", "TD", "traits", "design",
                                "spatial", "engine", "predicted", "sumTab"))
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
  modelSp <- fitTD(testTD, trials = "E1", design = "rowcol", traits = "t1")
  expect_SSA(modelSp)
  expect_SSAMod(modelSp, "mRand")
  expect_SSAMod(modelSp, "mFix")
  expect_identical(modelSp[["E1"]]$traits, "t1")
  expect_identical(modelSp[["E1"]]$design, "rowcol")
  expect_identical(modelSp[["E1"]]$spatial[["t1"]], "2 dimensional P-splines")
  expect_identical(modelSp[["E1"]]$engine, "SpATS")
})

test_that("running models creates objects with correct structure - lme4", {
  modelLm <- fitTD(testTD, trials = "E1", design = "rcbd", traits = "t1",
                   engine = "lme4")
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
  modelAs <- fitTD(testTD, trials = "E1", design = "res.ibd",
                   traits = "t1", engine = "asreml")
  expect_SSA(modelAs)
  expect_SSAMod(modelAs, "mRand")
  expect_SSAMod(modelAs, "mFix")
  expect_identical(modelAs[["E1"]]$traits, "t1")
  expect_identical(modelAs[["E1"]]$design, "res.ibd")
  expect_false(modelAs[["E1"]]$spatial[["t1"]])
  expect_identical(modelAs[["E1"]]$engine, "asreml")
})

test_that("option what produces expected output - SpATS", {
  modelSp <- fitTD(testTD, trials = "E1", design = "res.ibd", traits = "t1")
  modelSpF <- fitTD(testTD, trials = "E1", design = "res.ibd",
                    traits = "t1", what = "fixed")
  expect_SSA(modelSpF)
  expect_null(modelSpF[["E1"]]$mRand)
  expect_SSAMod(modelSpF, "mFix")
  expect_equal(modelSpF[["E1"]]$mFix, modelSp[["E1"]]$mFix)
  modelSpR <- fitTD(testTD, trials = "E1", design = "res.ibd",
                    traits = "t1", what = "random")
  expect_SSA(modelSpR)
  expect_SSAMod(modelSpR, "mRand")
  expect_equal(modelSp[["E1"]]$mRand, modelSpR[["E1"]]$mRand)
  expect_null(modelSpR[["E1"]]$mFix)
})

test_that("option what produces expected output - lme4", {
  modelLm <- fitTD(testTD, trials = "E1", design = "rcbd", traits = "t1",
                   engine = "lme4")
  modelLmF <- fitTD(testTD, trials = "E1", design = "rcbd",
                    traits = "t1", what = "fixed", engine = "lme4")
  expect_SSA(modelLmF)
  expect_null(modelLmF[["E1"]]$mRand)
  expect_SSAMod(modelLmF, "mFix", "lm")
  expect_equal(modelLmF[["E1"]]$mFix, modelLm[["E1"]]$mFix)
  modelLmR <- fitTD(testTD, trials = "E1", design = "rcbd",
                    traits = "t1", what = "random", engine = "lme4")
  expect_SSA(modelLmR)
  expect_SSAMod(modelLmR, "mRand", "lmerMod")
  expect_equal(modelLm[["E1"]]$mRand, modelLmR[["E1"]]$mRand)
  expect_null(modelLmR[["E1"]]$mFix)
})

test_that("option what produces expected output - asreml", {
  skip_on_cran()
  modelAs <- fitTD(testTD, trials = "E1", design = "rowcol", traits = "t1",
                   engine = "asreml")
  modelAsF <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = "t1", what = "fixed", engine = "asreml")
  expect_SSA(modelAsF)
  expect_null(modelAsF[["E1"]]$mRand)
  expect_SSAMod(modelAsF, "mFix")
  modelAsR <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = "t1", what = "random", engine = "asreml")
  expect_SSA(modelAsR)
  expect_SSAMod(modelAsR, "mRand")
  expect_null(modelAsR[["E1"]]$mFix)
})

test_that("running models for multiple traits produces correct output structure", {
  modelSp2 <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = paste0("t", 1:2))
  modelSp3 <- fitTD(testTD, trials = "E1", design = "rcbd",
                    traits = paste0("t", 1:3))
  modelSp4 <- fitTD(testTD, trials = "E1", design = "rowcol",
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
  modelSp1 <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = "t1")
  modelSp2 <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = paste0("t", 1:2))
  modelSp3 <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = paste0("t", 1:3))
  modelSp4 <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = paste0("t", 1:4))
  expect_equal(modelSp1[["E1"]]$mRand, modelSp2[["E1"]]$mRand["t1"])
  expect_equal(modelSp2[["E1"]]$mRand, modelSp3[["E1"]]$mRand[paste0("t", 1:2)])
  expect_equal(modelSp3[["E1"]]$mRand, modelSp4[["E1"]]$mRand[paste0("t", 1:3)])
  expect_equal(modelSp1[["E1"]]$mFix, modelSp2[["E1"]]$mFix["t1"])
  expect_equal(modelSp2[["E1"]]$mFix, modelSp3[["E1"]]$mFix[paste0("t", 1:2)])
  expect_equal(modelSp3[["E1"]]$mFix, modelSp4[["E1"]]$mFix[paste0("t", 1:3)])
})

test_that("option covariates produces expected output structure", {
  modelSpCov <- fitTD(testTD, trials = "E1", design = "rowcol",
                      traits = "t1", covariates = "repId")
  expect_SSA(modelSpCov)
  expect_SSAMod(modelSpCov, "mRand")
  expect_SSAMod(modelSpCov, "mFix")
  expect_true(grepl(pattern = "repId",
                    x = deparse(modelSpCov[["E1"]]$mRand$t1$model$fixed)))
  expect_true(grepl(pattern = "repId",
                    x = deparse(modelSpCov[["E1"]]$mFix$t1$model$fixed)))
  modelLmCov <- fitTD(testTD, trials = "E1", design = "rcbd",
                      traits = "t1", covariates = "repId", engine = "lme4")
  expect_SSA(modelLmCov)
  expect_SSAMod(modelLmCov, "mRand", "lmerMod")
  expect_SSAMod(modelLmCov, "mFix", "lm")
  expect_true("repId" %in% colnames(modelLmCov[["E1"]]$mRand$t1@frame))
  expect_true("repId" %in% colnames(modelLmCov[["E1"]]$mFix$t1$model))
  skip_on_cran()
  modelAsCov <- fitTD(testTD, trials = "E1", design = "rowcol",
                      traits = "t1", covariates = "repId",
                      engine = "asreml")
  expect_SSA(modelAsCov)
  expect_SSAMod(modelAsCov, "mRand")
  expect_SSAMod(modelAsCov, "mFix")
  if (asreml4()) {
    expect_true(grepl(pattern = "repId",
                      x = deparse(modelAsCov[["E1"]]$mRand$t1$formulae$fixed)))
    expect_true(grepl(pattern = "repId",
                      x = deparse(modelAsCov[["E1"]]$mFix$t1$formulae$fixed)))
  } else {
    expect_true(grepl(pattern = "repId",
                      x = deparse(modelAsCov[["E1"]]$mRand$t1$fixed.formula)))
    expect_true(grepl(pattern = "repId",
                      x = deparse(modelAsCov[["E1"]]$mFix$t1$fixed.formula)))
  }
})

test_that("option useCheckId produces expected output structure", {
  modelSpCi <- fitTD(testTD, trials = "E1", design = "rowcol", traits = "t1",
                     useCheckId = TRUE)
  expect_SSA(modelSpCi)
  expect_SSAMod(modelSpCi, "mRand")
  expect_SSAMod(modelSpCi, "mFix")
  expect_true(grepl(pattern = "checkId",
                    x = deparse(modelSpCi[["E1"]]$mRand$t1$model$fixed)))
  expect_true(grepl(pattern = "checkId",
                    x = deparse(modelSpCi[["E1"]]$mFix$t1$model$fixed)))
  modelLmCi <- fitTD(testTD, trials = "E1", design = "rcbd", traits = "t1",
                     useCheckId = TRUE, engine = "lme4")
  expect_SSA(modelLmCi)
  expect_SSAMod(modelLmCi, "mRand", "lmerMod")
  expect_SSAMod(modelLmCi, "mFix", "lm")
  expect_true("checkId" %in% colnames(modelLmCi[["E1"]]$mRand$t1@frame))
  expect_true("checkId" %in% colnames(modelLmCi[["E1"]]$mFix$t1$model))
})

test_that("option trySpatial produces expected output structure", {
  modelSp <- fitTD(testTD, trials = "E1", design = "rowcol", traits = "t1")
  modelSpTs <- fitTD(testTD, trials = "E1", design = "rowcol",
                     traits = "t1", trySpatial = TRUE)
  expect_SSA(modelSpTs)
  expect_SSAMod(modelSpTs, "mRand")
  expect_SSAMod(modelSpTs, "mFix")
  ## SpATS should use trySpatial as default. Timestamp will be different.
  expect_equivalent(modelSp, modelSpTs)
  expect_warning(fitTD(testTD, trials = "E1", design = "rowcol",
                       traits = "t1", trySpatial = TRUE, engine = "lme4"),
                 "Spatial models can only be fitted using SpATS or asreml.")
  skip_on_cran()
  modelAsTs <- fitTD(testTD, trials = "E1", design = "ibd", traits = "t1",
                     trySpatial = TRUE, engine = "asreml")
  expect_SSA(modelAsTs)
  expect_SSAMod(modelAsTs, "mRand")
  expect_SSAMod(modelAsTs, "mFix")
  expect_is(modelAsTs[["E1"]]$spatial, "list")
  expect_length(modelAsTs[["E1"]]$spatial, 1)
  expect_named(modelAsTs[["E1"]]$spatial, "t1")
  expect_identical(modelAsTs[["E1"]]$spatial$t1, "none")
})

test_that("option nSeg in control produces correct output", {
  ## Test using equivalence because of timestamp.
  modelSp <- fitTD(testTD, trials = "E1", design = "rowcol", traits = "t1",
                   control = list(nSeg = 1))
  modelSp1 <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = "t1", control = list(nSeg = c(1, 1)))
  expect_equivalent(modelSp, modelSp1)
  expect_error(fitTD(testTD, trials = "E1", design = "rowcol",
                     traits = "t1", control = list(nSeg = list(c(1, 1)))),
               "should be a named item in list of nSeg")
  modelSp2 <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = "t1",
                    control = list(nSeg = list(E1 = c(1, 1))))
  expect_equivalent(modelSp, modelSp2)
  modelSp3 <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = "t1",
                    control = list(nSeg = list(E3 = c(1, 1),
                                               E1 = c(1, 1))))
  expect_equivalent(modelSp, modelSp3)
})

test_that("option nestDiv in control produces correct output", {
  ## Test using equivalence because of timestamp.
  modelSp <- fitTD(testTD, trials = "E1", design = "rowcol", traits = "t1",
                   control = list(nestDiv = 3))
  modelSp1 <- fitTD(testTD, trials = "E1", design = "rowcol",
                    traits = "t1", control = list(nestDiv = c(3, 3)))
  expect_equivalent(modelSp, modelSp1)
  expect_warning(fitTD(testTD, trials = "E1", design = "rowcol",
                       traits = "t1", control = list(nestDiv = 0)),
                 "Invalid value for control parameter nestDiv")
})

test_that("option progress functions properly", {
  expect_output(fitTD(testTD, trials = "E1", design = "rowcol", traits = "t1",
                      progress = TRUE), "Fitting models for t1 in E1")
})

testData2 <- testData
## Set all observations to NA for 1 trial in 1 field to create data
## that causes the model engines to crash.
## STTunModel should be able to handle this and still produce output for
## the other models.
testData2[testData2$field == "E1", "t2"] <- NA
testTD2 <- createTD(data = testData2, trial = "field",
                    genotype = "seed", rowCoord = "Y", colCoord = "X")
test_that("Trial with missing data is handled properly when fitting models", {
  expect_warning(modelSp <- fitTD(testTD2, trials = "E1",
                                  design = "rowcol",
                                  traits = c("t1", "t2", "t3")),
                 "Error in SpATS")
  expect_SSA(modelSp)
  expect_warning(modelLm <- fitTD(testTD2, trials = "E1", design = "rowcol",
                                  traits = c("t1", "t2", "t3"),
                                  engine = "lme4"),
                 "Error in lmer")
  expect_SSA(modelLm)
  skip_on_cran()
  expect_warning(modelAs <- fitTD(testTD2, trials = "E1", design = "rowcol",
                                  traits = c("t1", "t2", "t3"),
                                  engine = "asreml"),
                 "Error in asreml")
  expect_SSA(modelAs)
  expect_warning(modelAs2 <- fitTD(testTD2, trials = "E1", design = "rowcol",
                                   traits = c("t1", "t2", "t3"),
                                   engine = "asreml", trySpatial = TRUE),
                 "Error in asreml")
  expect_SSA(modelAs2)
})

testData3 <- testData
## Set replicates to 1 for 1 field to test that design is changed to
## corresponding design without replicates.
testData3[testData3$field == "E1", "rep"] <- 1
testTD3 <- createTD(data = testData3, trial = "field",
                    genotype = "seed", rowCoord = "Y", colCoord = "X",
                    repId = "rep")
test_that("Design is modified when replicates contain only 1 distinct value", {
  expect_warning(modelSp <- fitTD(testTD3, trials = "E1",
                                  design = "res.rowcol", traits = c("t1")),
                 "Design changed")
  expect_equal(modelSp$E1$design, "rowcol")
})
