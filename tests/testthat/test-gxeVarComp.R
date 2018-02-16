context("gxeVarComp")

testTD <- createTD(data = testData, genotype = "seed",
                   env = "field", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoordinates = "Y", colCoordinates = "X")
BLUEsList <- lapply(X = levels(testTD$env), FUN = function(e) {
  modelSp <- STRunModel(testTD[testTD$env == e, ], design = "rowcol",
                        traits = c("t1", "t2", "t3", "t4"))
  STExtract(modelSp, what = "BLUEs", keep = "env")
})
BLUEs <- createTD(Reduce(f = rbind, x = BLUEsList))

geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
geVCLm <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4")
for (geVC in list(geVCAs, geVCLm)) {
  test_that("output is of the right class", {
    expect_is(geVC, "varComp")
    expect_is(geVC$model, "SSA")
    expect_is(geVC$choice, "character")
    expect_is(geVC$summary, "matrix")
    expect_is(geVC$vcov, "matrix")
    expect_is(geVC$criterion, "character")
    expect_is(geVC$engine, "character")
  })
}

summAs <- geVCAs$summary
test_that("asreml model gives correct output", {
  expect_equal(geVCAs$choice, "identity")
  expect_equal(rownames(summAs),
               c("identity", "cs", "diagonal", "outside", "hcs", "unstructured",
                 "fa", "fa2"))
  expect_equivalent(summAs[, "AIC"],
                    c(316.133524782187, 317.509091360267, 319.974781981561,
                      321.37350203058, 321.400027844011, 324.919291077454, NA, NA))
  expect_equivalent(summAs[, "BIC"],
                    c(317.534722163849, 320.311486123591, 324.178374126547,
                      326.978291557229, 327.00481737066, 333.326475367427, NA, NA))
  expect_equivalent(summAs[, "Deviance"],
                    c(314.133524782187, 313.509091360267, 313.974781981561,
                      313.37350203058, 313.400027844011, 312.919291077454, NA, NA))
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
  expect_equivalent(geVCLm$vcov, c(35.7978301877475, 4.51026038047197, 4.51026038047198,
                                   4.51026038047197, 35.7978301877475, 4.51026038047197,
                                   4.51026038047198, 4.51026038047197, 35.7978301877475))
})

geVCAsA <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml",
                      criterion = "AIC")
test_that("option criterion works properly", {
  expect_identical(geVCAs$summary[, "BIC"],
                   geVCAs$summary[, "BIC"][order(geVCAs$summary[, "BIC"])])
  expect_identical(geVCAsA$summary[, "AIC"],
                   geVCAsA$summary[, "AIC"][order(geVCAsA$summary[, "AIC"])])
})
