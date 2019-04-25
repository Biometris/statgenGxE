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
  expect_equivalent(sumSp$heritability, 0.46)
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
  expect_equivalent(sumLm$heritability, 0.121944484712314)
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
  expect_equivalent(sumAs$heritability, 0.206084437486658)
  expect_equal(nrow(sumAs$sed), 3)
  expect_equal(nrow(sumAs$lsd), 3)
})

test_that("summary.SSA produces correct output for multiple trials", {
  modelSp <- STRunModel(testTD, trials = c("E1", "E2"), design = "rowcol",
                        traits = "t1")
  sumSp <- summary(modelSp, traits = "t1")
  expect_length(sumSp, 3)
  expect_is(sumSp$sumTab, "matrix")
  expect_equal(dim(sumSp$sumTab), c(2, 9))
  expect_equal(sumSp$what, "BLUEs")
})

test_that("option nBest functions properly", {
  sumSp <- summary(modelSp, nBest = 5)
  expect_equal(dim(sumSp$meanTab), c(5, 4))
})

test_that("print.summary.SSA functions properly", {
  sumSp <- capture.output(print(summary(modelSp)))
  sumSp2 <- capture.output(print(summary(modelSp, nBest = NA)))
  expect_true(all(c("Summary statistics for t1 in E1  ",
                    "Estimated heritability ",
                    "Predicted means (BLUEs & BLUPs) ") %in% sumSp))
  expect_false(any(grepl("Best", sumSp2)))
  skip_on_cran()
  modelAs <- STRunModel(testTD, trials = "E1", design = "rowcol", traits = "t1",
                        engine = "asreml")
  sumAs <- capture.output(print(summary(modelAs)))
  expect_true(all(c("Standard Error of Difference (genotype modeled as fixed effect) ",
                    "Least Significant Difference (genotype modeled as fixed effect) ") %in%
                    sumAs))
})

test_that("print.summary.SSA functions properly for multiple trials", {
  modelSp <- STRunModel(testTD, trials = c("E1", "E2"), design = "rowcol",
                        traits = "t1")
  sumSp <- capture.output(print(summary(modelSp)))
  expect_true("Summary statistics for BLUEs of t1 " %in% sumSp)
})

test_that("function SSAtoTD functions properly", {
  TDSp <- SSAtoTD(SSA = modelSp)
  expect_is(TDSp, "TD")
  expect_equal(colnames(TDSp$E1),
               c("genotype", "trial", "BLUEs_t1", "seBLUEs_t1", "BLUPs_t1",
                 "seBLUPs_t1"))
  TDSp2 <- SSAtoTD(SSA = modelSp, what = "BLUEs")
  expect_equal(colnames(TDSp2$E1),
               c("genotype", "trial", "t1"))
  expect_warning(TDSp3 <- SSAtoTD(SSA = modelSp, what = "BLUEs", addWt = TRUE),
                 "Weights can only be added together with seBLUEs")
  expect_equal(colnames(TDSp3$E1),
               c("genotype", "trial", "BLUEs_t1", "seBLUEs_t1", "wt"))
  expect_warning(TDSp4 <- SSAtoTD(SSA = modelSp, keep = "family"),
                 "Duplicate values for")
  expect_equal(colnames(TDSp4$E1),
               c("genotype", "trial", "BLUEs_t1", "seBLUEs_t1", "BLUPs_t1",
                 "seBLUPs_t1"))
})

test_that("function SSAtoCross functions properly", {
  myModel <- STRunModel(TD = TDHeat05, design = "res.rowcol", traits = "yield",
                        what = "fixed")
  expect_error(SSAtoCross(SSA = myModel, trial = "HEAT06",
                          genoFile = system.file("extdata", "markers.csv",
                                                 package = "RAP")),
               "single character string defining a trial in SSA")
  cross <- SSAtoCross(SSA = myModel,
                      genoFile = system.file("extdata", "markers.csv",
                                             package = "RAP"))
  expect_is(cross, "cross")
  expect_is(cross$pheno, "data.frame")
  expect_equal(dim(cross$pheno), c(169, 2))
})

test_that("function report.SSA functions properly" ,{
  ## Reporting doesn't work on cran because of usage of pdflatex.
  skip_on_cran()
  modelSp <- STRunModel(testTD, trials = c("E1", "E2"), design = "rowcol",
                        traits = c("t1", "t2"))
  expect_error(report(modelSp),
               "No trial provided but multiple trials found")
  expect_error(report(modelSp, trial = "E3"),
               "single character string defining a trial in SSA")
  expect_error(report(modelSp, trial = "E1"),
               "No trait provided but multiple traits found")
  expect_error(report(modelSp, trial = "E1", trait = "t3"),
               "single character string defining a trait for which a model")
  tmpFile = tempfile(fileext = ".pdf")
  expect_silent(report(modelSp, trial = "E1", trait = "t1", outfile = tmpFile))
  expect_true(file.exists(tmpFile))
  expect_true(file.exists(gsub(pattern = ".pdf", replacement = ".tex",
                               x = tmpFile)))
})
