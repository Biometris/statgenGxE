context("gxeMegaEnv")

testTD <- createTD(data = testData, genotype = "seed",
                   env = "field", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoordinates = "Y", colCoordinates = "X")

BLUEsList <- lapply(X = levels(testTD$env), FUN = function(e) {
  modelSp <- STRunModel(testTD[testTD$env == e, ], design = "rowcol",
                        traits = "t1")
  STExtract(modelSp, what = "BLUEs", keep = "env")
})
BLUEs <- createTD(Reduce(f = rbind, x = BLUEsList))

geMegaEnv <- gxeMegaEnv(TD = BLUEs, trait = "t1")
test_that("mega-environments are computed correctly", {
  expect_is(geMegaEnv, "TD")
  expect_identical(dim(geMegaEnv), c(45L, 4L))
  expect_equal(as.numeric(geMegaEnv$megaEnv), rep(x = c(1, 2, 3), each = 15))
  expect_equal(levels(geMegaEnv$megaEnv), c("3", "2", "1"))
})

summ <- attr(x = geMegaEnv, which = "sumTab")
test_that("summary is computed correctly", {
  expect_is(summ, "data.frame")
  expect_equal(as.character(summ$winGeno), c("G14", "G4", "G5"))
  expect_equal(summ$`AMMI estimates`,
               c(137.826424971715, 109.450972152511, 127.386419996042))
  expect_null(attr(x = gxeMegaEnv(TD = BLUEs, trait = "t1", sumTab = FALSE),
              which = "sumTab"))
})

geMegaEnvMin <- gxeMegaEnv(TD = BLUEs, trait =
                             "t1", method = "min")
test_that("option method functions properly", {
  expect_equal(as.numeric(geMegaEnv$megaEnv), rep(x = c(1, 2, 3), each = 15))
  expect_equal(levels(geMegaEnv$megaEnv), c("3", "2", "1"))
  expect_equal(as.character(attr(x = geMegaEnvMin, which = "sumTab")$winGeno),
              c("G12", "G15", "G6"))
  expect_equal(attr(x = geMegaEnvMin, which = "sumTab")$`AMMI estimates`,
               c(55.0432001091191, 40.7120064386264, 35.3752613055704))
})
