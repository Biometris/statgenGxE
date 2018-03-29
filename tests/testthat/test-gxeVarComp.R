context("gxeVarComp")

testTD <- createTD(data = testData, genotype = "seed",
                   trial = "field", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoordinates = "Y", colCoordinates = "X")
modelSp <- STRunModel(testTD, design = "rowcol",
                      traits = c("t1", "t2", "t3", "t4"))
BLUEsList <- STExtract(modelSp, what = "BLUEs", keep = "trial")
BLUEs <- createTD(Reduce(f = rbind, x = BLUEsList))

geVCLm <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4")
test_that("output is of the right class for lme4", {
  expect_is(geVCLm, "varComp")
  expect_is(geVCLm$SSA, "SSA")
  expect_is(geVCLm$choice, "character")
  expect_is(geVCLm$summary, "matrix")
  expect_is(geVCLm$vcov, "matrix")
  expect_is(geVCLm$criterion, "character")
  expect_is(geVCLm$engine, "character")
})

if (requireNamespace("asreml", quietly = TRUE)) {
  geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
}
test_that("output is of the right class for asreml", {
  testthat::skip_on_cran()
  expect_is(geVCAs, "varComp")
  expect_is(geVCAs$SSA, "SSA")
  expect_is(geVCAs$choice, "character")
  expect_is(geVCAs$summary, "matrix")
  expect_is(geVCAs$vcov, "matrix")
  expect_is(geVCAs$criterion, "character")
  expect_is(geVCAs$engine, "character")
})

test_that("asreml model gives correct output", {
  testthat::skip_on_cran()
  summAs <- geVCAs$summary
  expect_equal(geVCAs$choice, "identity")
  expect_equal(rownames(summAs),
               c("identity", "cs", "diagonal", "outside", "hcs", "unstructured",
                 "fa", "fa2"))
  expect_equivalent(summAs[, "AIC"],
                    c(309.722795184198, 311.145684341532, 313.533669549437,
                      314.991580913606, 315.030761318384, 318.54873092037, NA, NA))
  expect_equivalent(summAs[, "BIC"],
                    c(311.12399256586, 313.948079104857, 317.737261694423,
                      320.596370440255, 320.635550845033, 326.955915210343, NA, NA))
  expect_equivalent(summAs[, "Deviance"],
                    c(307.722795184198, 307.145684341532, 307.533669549437,
                      306.991580913606, 307.030761318384, 306.54873092037, NA, NA))
  expect_equivalent(summAs[, "NParameters"], c(1, 2, 3, 4, 4, 6, NA, NA))
  expect_equivalent(geVCAs$vcov, c(36.6317945197525, 0, 0, 0, 36.6317945197525,
                                   0, 0, 0, 39.2483512711634))
})

summLm <- geVCLm$summary
test_that("lme4 model gives correct output", {
  expect_equal(geVCLm$choice, "cs")
  expect_equal(rownames(summLm), "cs")
  expect_equivalent(summLm[, "AIC"], 386.498644064315)
  expect_equivalent(summLm[, "BIC"], 390.067023332152)
  expect_equivalent(summLm[, "Deviance"], 382.498644064315)
  expect_equivalent(summLm[, "NParameters"], 2)
  ## This test works fine in RStudio but gives an error when testing on CRAN.
  ## Therefore added a lower tolerance
  expect_equivalent(geVCLm$vcov,
                    c(36.6109167394648, 4.44290083464202, 4.442900834642,
                      4.44290083464202, 36.6109167394648, 4.442900834642,
                      4.442900834642, 4.442900834642, 39.1572940974083),
                    tolerance = 1e-6)
})

test_that("option criterion works properly", {
  testthat::skip_on_cran()
  geVCAsA <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml",
                        criterion = "AIC")
  expect_identical(geVCAs$summary[, "BIC"],
                   geVCAs$summary[, "BIC"][order(geVCAs$summary[, "BIC"])])
  expect_identical(geVCAsA$summary[, "AIC"],
                   geVCAsA$summary[, "AIC"][order(geVCAsA$summary[, "AIC"])])
})
