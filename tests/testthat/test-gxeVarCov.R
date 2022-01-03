context("gxeVarCov")

test_that("general checks in gxeVarCov function properly", {
  expect_error(gxeVarCov(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeVarCov(BLUEs, trait = "t5"),
               "t5 has to be a column in TD")
  expect_error(gxeVarCov(BLUEs, trials = "E4", trait = "t1"),
               "a character vector defining trials in BLUEs")
})

geVCLm <- gxeVarCov(TD = BLUEs, trait = "t1", engine = "lme4")
test_that("output is of the right class for lme4", {
  expect_is(geVCLm, "varCov")
  expect_is(geVCLm$STA, "STA")
  expect_is(geVCLm$choice, "character")
  expect_is(geVCLm$summary, "matrix")
  expect_is(geVCLm$vcov, "matrix")
  expect_is(geVCLm$criterion, "character")
  expect_is(geVCLm$engine, "character")
})

test_that("output is of the right class for asreml", {
  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarCov(TD = BLUEs, trait = "t1", engine = "asreml")
  expect_is(geVCAs, "varCov")
  expect_is(geVCAs$STA, "STA")
  expect_is(geVCAs$choice, "character")
  expect_is(geVCAs$summary, "matrix")
  expect_is(geVCAs$vcov, "matrix")
  expect_is(geVCAs$criterion, "character")
  expect_is(geVCAs$engine, "character")
})

test_that("asreml model gives correct output", {
  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarCov(TD = BLUEs, trait = "t1", engine = "asreml")
  summAs <- geVCAs$summary
  expect_equal(geVCAs$choice, "identity")
  expect_equal(rownames(summAs),
               c("identity", "cs", "diagonal", "hcs", "outside", "unstructured",
                 "fa", "fa2"))
  expect_equivalent(summAs[, "AIC"],
                    c(302.538128471826, 302.958375388098, 306.507199541997,
                      306.866481338154, 308.538032578925, 308.417501699175,
                      NA, NA))
  expect_equivalent(summAs[, "BIC"],
                    c(303.939325853488, 305.760770151422, 310.710791686983,
                      312.471270864803, 314.142822105573, 316.824685989148,
                      NA, NA))
  expect_equivalent(summAs[, "Deviance"],
                    c(300.538128471826, 298.958375388098, 300.507199541997,
                      298.866481338154, 300.538032578925, 296.417501699175,
                      NA, NA))
  expect_equivalent(summAs[, "NParameters"], c(1, 2, 3, 4, 4, 6, NA, NA))
  expect_equivalent(geVCAs$vcov, c(25.8985599609355, 0, 0, 0, 25.8985599609355,
                                   0, 0, 0, 25.8985599609355))
})

test_that("lme4 model gives correct output", {
  summLm <- geVCLm$summary
  expect_equal(geVCLm$choice, "cs")
  expect_equal(rownames(summLm), "cs")
  expect_equivalent(summLm[, "AIC"], 380.14921134723)
  expect_equivalent(summLm[, "BIC"], 383.762536326771)
  expect_equivalent(summLm[, "Deviance"], 376.14921134723)
  expect_equivalent(summLm[, "NParameters"], 2)
  expect_equivalent(geVCLm$vcov,
                    c(25.8985599502474, 5.25233386534017, 5.25233386534018,
                      5.25233386534017, 25.8985599502474, 5.25233386534018,
                      5.25233386534018, 5.25233386534018, 25.8985599502474),
                    tolerance = 1e-6)
})

test_that("option criterion works properly", {
  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarCov(TD = BLUEs, trait = "t1", engine = "asreml")
  geVCAsA <- gxeVarCov(TD = BLUEs, trait = "t1", engine = "asreml",
                        criterion = "AIC")
  expect_equal(geVCAs$summary[, "BIC"],
               geVCAs$summary[, "BIC"][order(geVCAs$summary[, "BIC"])])
  expect_equal(geVCAsA$summary[, "AIC"],
               geVCAsA$summary[, "AIC"][order(geVCAsA$summary[, "AIC"])])
})

test_that("models for fa and fa2 are fitted when #trials >= 5", {
  skip_on_cran()
  skip_on_ci()
  expect_warning(capture_output(geVC <- gxeVarCov(TD = BLUEsYear, trait = "t1",
                                                   engine = "asreml")))
  expect_equal(rownames(geVC$summary),
               c("identity", "cs", "diagonal", "hcs", "fa", "fa2",
                 "unstructured", "outside"))
})

test_that("option models works properly", {
  skip_on_cran()
  skip_on_ci()
  expect_error(gxeVarCov(TD = BLUEs, trait = "t1", models = "fa"),
               "With lme4 only compound symmetry")
  geVCMod <- gxeVarCov(TD = BLUEs, trait = "t1",
                      models = c("identity", "cs", "fa"),
                      engine = "asreml")
  summMod <- geVCMod$summary
  expect_equal(rownames(summMod), c("identity", "cs", "fa"))
})
