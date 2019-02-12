context("gxeVarComp")

testTD <- createTD(data = testData, genotype = "seed",
                   trial = "field", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoord = "Y", colCoord = "X")
modelSp <- STRunModel(testTD, design = "rowcol",
                      traits = c("t1", "t2", "t3", "t4"))
BLUEs <- SSAtoTD(modelSp, what = "BLUEs")

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
               c("identity", "cs", "diagonal", "hcs", "outside", "unstructured",
                 "fa", "fa2"))
  expect_equivalent(summAs[, "AIC"],
                    c(316.133524782187, 317.509091360267, 319.974781981561,
                      321.400027844011, 322.133464280762, 324.919289306661,
                      NA, NA))
  expect_equivalent(summAs[, "BIC"],
                    c(317.534722163849, 320.311486123591, 324.178374126547,
                      327.00481737066, 327.73825380741, 333.326473596634,
                      NA, NA))
  expect_equivalent(summAs[, "Deviance"],
                    c(314.133524782187, 313.509091360267, 313.974781981561,
                      313.400027844011, 314.133464280762, 312.919289306661,
                      NA, NA))
  expect_equivalent(summAs[, "NParameters"], c(1, 2, 3, 4, 4, 6, NA, NA))
  expect_equivalent(geVCAs$vcov, c(35.7978300824404, 0, 0, 0, 35.7978300824404,
                                   0, 0, 0, 35.7978300824404))
})

summLm <- geVCLm$summary
test_that("lme4 model gives correct output", {
  expect_equal(geVCLm$choice, "cs")
  expect_equal(rownames(summLm), "cs")
  expect_equivalent(summLm[, "AIC"], 394.699928149459)
  expect_equivalent(summLm[, "BIC"], 398.313253129)
  expect_equivalent(summLm[, "Deviance"], 390.699928149459)
  expect_equivalent(summLm[, "NParameters"], 2)
  ## This test works fine in RStudio but gives an error when testing on CRAN.
  ## Therefore added a lower tolerance
  expect_equivalent(geVCLm$vcov,
                    c(35.7978169002067, 4.5102010045834, 4.5102010045834,
                      4.5102010045834, 35.7978169002067, 4.5102010045834,
                      4.5102010045834, 4.5102010045834, 35.7978169002067),
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
