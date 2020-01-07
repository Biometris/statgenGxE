context("gxeVarComp")

modelSp <- fitTD(testTD, design = "rowcol", traits = "t1")
BLUEs <- STAtoTD(modelSp, what = "BLUEs")

modelSpYear <- fitTD(testTDYear, design = "rowcol", traits = "t1")
BLUEsYear <- STAtoTD(modelSpYear, what = "BLUEs")

test_that("general checks in gxeAmmi function properly", {
  expect_error(gxeVarComp(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeVarComp(BLUEs, trait = "t5"),
               "t5 has to be a column in TD")
  expect_error(gxeVarComp(BLUEs, trials = "E4", trait = "t1"),
               "a character vector defining trials in BLUEs")
})

geVCLm <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4")
test_that("output is of the right class for lme4", {
  expect_is(geVCLm, "varComp")
  expect_is(geVCLm$STA, "STA")
  expect_is(geVCLm$choice, "character")
  expect_is(geVCLm$summary, "matrix")
  expect_is(geVCLm$vcov, "matrix")
  expect_is(geVCLm$criterion, "character")
  expect_is(geVCLm$engine, "character")
})

test_that("varComp models are fitted correctly", {
  skip_on_cran()
  expect_warning(gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml"),
                 "No convergence for outside")
  expect_warning(capture_output(gxeVarComp(TD = BLUEsYear, trait = "t1",
                                           engine = "asreml")),
                 "Asreml gave the following error for outside")
})

test_that("output is of the right class for asreml", {
  skip_on_cran()
  expect_warning(geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1",
                                      engine = "asreml"))
  expect_is(geVCAs, "varComp")
  expect_is(geVCAs$STA, "STA")
  expect_is(geVCAs$choice, "character")
  expect_is(geVCAs$summary, "matrix")
  expect_is(geVCAs$vcov, "matrix")
  expect_is(geVCAs$criterion, "character")
  expect_is(geVCAs$engine, "character")
})

test_that("asreml model gives correct output", {
  skip_on_cran()
  expect_warning(geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1",
                                      engine = "asreml"))
  summAs <- geVCAs$summary
  expect_equal(geVCAs$choice, "identity")
  expect_equal(rownames(summAs),
               c("identity", "cs", "diagonal", "hcs", "unstructured", "outside",
                 "fa", "fa2"))
  expect_equivalent(summAs[, "AIC"],
                    c(316.133524782187, 317.509091360399, 319.974782019118,
                      321.400028597069, 324.919296781454, Inf, NA, NA))
  expect_equivalent(summAs[, "BIC"],
                    c(317.534722163849, 320.311486123591, 324.178374126547,
                      327.00481737066, 333.326473596634, Inf, NA, NA))
  expect_equivalent(summAs[, "Deviance"],
                    c(314.133524782187, 313.509091360267, 313.974781981561,
                      313.400027844011, 312.919289306661, Inf, NA, NA))
  expect_equivalent(summAs[, "NParameters"], c(1, 2, 3, 4, 6, 4, NA, NA))
  expect_equivalent(geVCAs$vcov, c(35.7978300824404, 0, 0, 0, 35.7978300824404,
                                   0, 0, 0, 35.7978300824404))
})

test_that("lme4 model gives correct output", {
  summLm <- geVCLm$summary
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
  skip_on_cran()
  expect_warning(geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1",
                                      engine = "asreml"))
  expect_warning(geVCAsA <- gxeVarComp(TD = BLUEs, trait = "t1",
                                       engine = "asreml",criterion = "AIC"))
  expect_identical(geVCAs$summary[, "BIC"],
                   geVCAs$summary[, "BIC"][order(geVCAs$summary[, "BIC"])])
  expect_identical(geVCAsA$summary[, "AIC"],
                   geVCAsA$summary[, "AIC"][order(geVCAsA$summary[, "AIC"])])
})

test_that("models for fa and fa2 are fitted when #trials >= 5", {
  skip_on_cran()
  expect_warning(capture_output(geVC <- gxeVarComp(TD = BLUEsYear, trait = "t1",
                                                   engine = "asreml")))
  expect_equal(rownames(geVC$summary),
               c("identity", "cs", "diagonal", "hcs", "fa", "fa2",
                 "unstructured", "outside"))
})
