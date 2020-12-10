context("gxeMegaEnv")

test_that("general checks in gxeMegaEnv function properly", {
  expect_error(gxeMegaEnv(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeMegaEnv(testTD, trait = c("t1", "t1")),
               "trait has to be a character string of length 1")
  expect_error(gxeMegaEnv(testTD, trait = "t5"),
               "t5 has to be a column in TD")
  expect_error(gxeMegaEnv(testTD, trait = "t1", byYear = TRUE),
               "year has to be a column in TD")
})

geMegaEnv <- gxeMegaEnv(TD = BLUEs, trait = "t1")
geMegaEnvTot <- Reduce(f = rbind, x = geMegaEnv$TD)
test_that("mega environments are computed correctly", {
  expect_is(geMegaEnv, "megaEnv")
  expect_is(geMegaEnv$TD, "TD")
  expect_is(geMegaEnv$summTab, "data.frame")
  expect_is(geMegaEnv$trait, "character")
  expect_equal(dim(geMegaEnvTot), c(45, 7))
  expect_equal(as.numeric(geMegaEnvTot[["megaEnv"]]),
               rep(x = c(2, 1, 1), each = 15))
  expect_equal(levels(geMegaEnvTot[["megaEnv"]]), c("megaEnv_1", "megaEnv_2"))
})

test_that("summary is computed correctly", {
  summ <- geMegaEnv$summTab
  expect_equal(as.character(summ[["Winning_genotype"]]), c("G2", "G2", "G9"))
  expect_equal(summ[["AMMI_estimates"]],
               c(111.48629093387, 119.206304548313, 126.941193543844))
})

geMegaEnvMin <- gxeMegaEnv(TD = BLUEs, trait = "t1", method = "min")
geMegaEnvMinTot <- Reduce(f = rbind, x = geMegaEnvMin$TD)
test_that("option method functions properly", {
  expect_equal(as.numeric(geMegaEnvMinTot[["megaEnv"]]),
               rep(x = c(3, 1, 2), each = 15))
  expect_equal(levels(geMegaEnvMinTot[["megaEnv"]]),
               c("megaEnv_1", "megaEnv_2", "megaEnv_3"))
  summMin <- geMegaEnvMin$summTab
  expect_equal(as.character(summMin[["Winning_genotype"]]),
               c("G10", "G6", "G7"))
  expect_equal(summMin[["AMMI_estimates"]],
               c(46.7063274979725, 53.6811189560119, 49.8026925088618))
})

test_that("existing megaEnv in TD is overwritten", {
  expect_warning(gxeMegaEnv(TD = geMegaEnv$TD, trait = "t1"),
                 "TD already contains a column megaEnv")
})

test_that("general checks in predict.megaEnv function properly", {
  expect_error(predict(geMegaEnv, trait = "t1", useYear = TRUE),
               "year has to be a column in TD")
})

test_that("predict.megaEnv functions correctly", {
  ## More random effects than observations, so empty data.frame returned.
  expect_warning(predict(geMegaEnvMin), "Empty data.frame returned")
  expect_warning(geTabLm <- predict(geMegaEnv),
                 "mega environments that are based on less than 10 trials")
  expect_is(geTabLm, "list")
  expect_length(geTabLm, 2)
  expect_named(geTabLm, c("predictedValue", "standardError"))
  expect_is(geTabLm$predictedValue, "data.frame")
  expect_is(geTabLm$standardError, "data.frame")
  expect_equivalent(geTabLm$predictedValue[1, ],
                    c(79.2407761061078, 79.4858110738655))
  expect_equivalent(geTabLm$standardError[1, ],
                    c(6.83087991095342, 6.38912904851271))
  skip_on_cran()
  skip_on_ci()
  expect_warning(geTabAs <- predict(geMegaEnv, engine = "asreml"),
                 "mega environments that are based on less than 10 trials")
  expect_equivalent(geTabAs$predictedValue[1, ],
                    c(79.2091290087036, 79.4263169941458))
  expect_equivalent(geTabAs$standardError[1, ],
                    c(7.20520150684489, 6.96252065379726))
})

## Modify data so it contains a year variable.
geMegaEnvNw <- list(TD = c(geMegaEnv$TD, geMegaEnv$TD), trait = "t1")
class(geMegaEnvNw) <- "megaEnv"

names(geMegaEnvNw$TD) <- paste0("E", 1:6)
class(geMegaEnvNw$TD) <- "TD"
geMegaEnvNw$TD$E1[["year"]] <- geMegaEnvNw$TD$E2[["year"]] <-
  geMegaEnvNw$TD$E3[["year"]] <- 1
geMegaEnvNw$TD$E4[["year"]] <- geMegaEnvNw$TD$E5[["year"]] <-
  geMegaEnvNw$TD$E6[["year"]] <- 2

test_that("option year in gxeTable functions properly", {
  expect_warning(megaEnvPred <- predict(geMegaEnvNw, useYear = TRUE),
                 "mega environments that are based on less than 10 trials")
  expect_equivalent(megaEnvPred$predictedValue[1, ],
                    c(75.7340426279225, 72.8307799806527))
  expect_equivalent(megaEnvPred$standardError[1, ],
                    c(5.69948967193805, 7.83055991157765))
})

test_that("option year in gxeTable functions properly for asreml", {
  skip_on_cran()
  skip_on_ci()
  expect_warning(megaEnvPredAs <- predict(geMegaEnvNw, useYear = TRUE,
                                          engine = "asreml"),
                 "mega environments that are based on less than 10 trials")
  expect_equivalent(megaEnvPredAs$predictedValue[1, ],
                    c(76.634397082771, 73.7311917123637))
  expect_equivalent(megaEnvPredAs$standardError[1, ],
                    c(5.87848501191572, 8.30017045501152))
})
