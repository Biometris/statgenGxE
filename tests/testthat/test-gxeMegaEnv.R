context("gxeMegaEnv")

testTD <- createTD(data = testData, genotype = "seed",
                   trial = "field", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoord = "Y", colCoord = "X")
modelSp <- STRunModel(testTD, design = "rowcol", traits = "t1")
BLUEs <- SSAtoTD(modelSp, what = "BLUEs")

tmp <- capture_output(geMegaEnv <- gxeMegaEnv(TD = BLUEs, trait = "t1"))
geMegaEnvTot <- Reduce(f = rbind, x = geMegaEnv)
test_that("mega-environments are computed correctly", {
  expect_is(geMegaEnv, "TD")
  expect_identical(dim(geMegaEnvTot), c(45L, 4L))
  expect_equal(as.numeric(geMegaEnvTot$megaEnv), rep(x = c(1, 2, 3), each = 15))
  expect_equal(levels(geMegaEnvTot$megaEnv), c("3", "2", "1"))
})

summ <- attr(x = geMegaEnv, which = "sumTab")
test_that("summary is computed correctly", {
  expect_is(summ, "data.frame")
  expect_equal(as.character(summ$`Winning genotype`), c("G14", "G4", "G5"))
  expect_equal(summ$`AMMI estimates`,
               c(137.826424971715, 109.450972152511, 127.386419996042))
})

tmp <- capture_output(geMegaEnvMin <- gxeMegaEnv(TD = BLUEs, trait = "t1",
                                                 method = "min"))
geMegaEnvMinTot <- Reduce(f = rbind, x = geMegaEnvMin)
test_that("option method functions properly", {
  expect_equal(as.numeric(geMegaEnvMinTot$megaEnv), rep(x = c(1, 2, 3),
                                                        each = 15))
  expect_equal(levels(geMegaEnvMinTot$megaEnv), c("3", "2", "1"))
  expect_equal(as.character(attr(x = geMegaEnvMin,
                                 which = "sumTab")$`Winning genotype`),
               c("G12", "G15", "G6"))
  expect_equal(attr(x = geMegaEnvMin, which = "sumTab")$`AMMI estimates`,
               c(55.0432001091191, 40.7120064386264, 35.3752613055704))
})

tmp <- capture_output(geMegaEnvMin <- gxeMegaEnv(TD = BLUEs, trait = "t1",
                                                 sumTab = FALSE))
test_that("option sumTab succesfully suppresses output", {
  expect_equal(tmp, "")
})

## Manipulate geMegaEnv to get proper output.
geMegaEnvNw <- geMegaEnv
geMegaEnvNw[["E3"]]$megaEnv <- 2
test_that("gxeTable functions correctly", {
  ## More random effects than observations, so empty data.frame returned.
  expect_warning(gxeTable(TD = geMegaEnv, trait = "t1"),
                 "Empty data.frame returned")
  geTabLm <- gxeTable(TD = geMegaEnvNw, trait = "t1")
  expect_equivalent(geTabLm$predictedValue[1, ],
                    c(80.4703105526859, 80.0424618298571))
  ## This test works fine in RStudio but gives an error when testing on CRAN.
  ## Therefore added a lower tolerance
  expect_equivalent(geTabLm$standardError[1, ],
                    c(9.50660071605727, 5.93819059539989), tolerance = 1e-7)
  testthat::skip_on_cran()
  expect_warning(gxeTable(TD = geMegaEnv, trait = "t1"),
                 "Empty data.frame returned")
  geTabAs <- gxeTable(TD = geMegaEnvNw, trait = "t1", engine = "asreml")
  expect_equivalent(geTabAs$predictedValue[1, ],
                    c(80.4703103942918, 80.0424619340866))
  expect_equivalent(geTabAs$standardError[1, ],
                    c(9.77290639163707, 6.58334962888759))
})



