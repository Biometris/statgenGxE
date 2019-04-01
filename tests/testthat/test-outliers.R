context("outliers")

## Create testdata containing only one trial.
testTD <- createTD(data = testData, trial = "field",
                   genotype = "seed", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoord = "Y", colCoord = "X")
modelLm <- STRunModel(testTD, trials = "E1", design = "rowcol",
                      traits = paste0("t", 1:4), engine = "lme4")

test_that("outliersSSA functions properly", {
  out1 <- outlierSSA(modelLm, trial = "E1", traits = "t1", verbose = FALSE)
  expect_is(out1, "list")
  expect_length(out1, 2)
  expect_is(out1$indicator, "data.frame")
  expect_equal(sum(out1$indicator[["t1"]]), 0)
  expect_null(out1$outliers)
})

test_that("outliersSSA functions properly for multiple traits", {
  out14 <- outlierSSA(modelLm, trial = "E1", traits = paste0("t", 1:4),
                      verbose = FALSE)
  expect_is(out14, "list")
  expect_length(out14, 2)
  expect_is(out14$indicator, "data.frame")
  expect_equal(colnames(out14$indicator), paste0("t", 1:4))
  expect_equivalent(colSums(out14$indicator), rep(x = 0, times = 4))
  expect_null(out14$outliers)
})

test_that("option what functions properly", {
  out1 <- outlierSSA(modelLm, trial = "E1", traits = "t1", what = "random",
                     verbose = FALSE)
  expect_is(out1, "list")
  expect_length(out1, 2)
  expect_is(out1$indicator, "data.frame")
  expect_equal(sum(out1$indicator[["t1"]]), 1)
  expect_equal(nrow(out1$outliers), 1)
  expect_equal(out1$outliers$res, 2.00884115446144)
})

test_that("option rLimit funtions properly", {
  out1 <- outlierSSA(modelLm, trial = "E1", traits = "t1", rLimit = 1,
                     verbose = FALSE)
  expect_equal(sum(out1$indicator[["t1"]]), 4)
  expect_equal(nrow(out1$outliers), 4)
  expect_equal(out1$outliers$res, c(-1.69959940078883, 1.16726345597009,
                                    1.69959940078883, -1.16726345597009))
})

test_that("option rLimit funtions properly for multiple traits", {
  out14 <- outlierSSA(modelLm, trial = "E1", traits = paste0("t", 1:4),
                      rLimit = 1, verbose = FALSE)
  expect_equivalent(colSums(out14$indicator), c(4, 6, 2, 4))
  expect_equal(nrow(out14$outliers), 16)
})

test_that("option commonFactors functions properly", {
  out1 <- outlierSSA(modelLm, trial = "E1", traits = "t1", rLimit = 1,
                     commonFactors = "subBlock", verbose = FALSE)
  expect_equal(sum(out1$indicator[["t1"]]), 4)
  expect_equal(nrow(out1$outliers), 12)
  expect_equal(sum(out1$outliers$similar), 8)
})

test_that("option verbose functions properly", {
  printOut1 <- capture.output(out1 <- outlierSSA(modelLm, trial = "E1",
                                                 traits = "t1", verbose = TRUE))
  printOut2 <- capture.output(out1 <- outlierSSA(modelLm, trial = "E1",
                                                 traits = "t1", what = "random",
                                                 verbose = TRUE))
  expect_equal(printOut1, "No large standardized residuals.")
  expect_true("Large standardized residuals" %in% printOut2)
  expect_true(any(grepl(pattern = "2.008841", x = printOut2)))
})

