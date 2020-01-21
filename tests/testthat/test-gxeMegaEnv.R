context("gxeMegaEnv")

modelSp <- fitTD(testTD, design = "rowcol", traits = "t1")
BLUEs <- STAtoTD(modelSp, what = "BLUEs")

test_that("general checks in gxeMegaEnv function properly", {
  expect_error(gxeMegaEnv(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeMegaEnv(testTD, trait = "t5"),
               "t5 has to be a column in TD")
  expect_error(gxeMegaEnv(testTD, trait = "t1", byYear = TRUE),
               "year has to be a column in TD")
})

geMegaEnv <- gxeMegaEnv(TD = BLUEs, trait = "t1", sumTab = FALSE)
geMegaEnvTot <- Reduce(f = rbind, x = geMegaEnv)
test_that("mega environments are computed correctly", {
  expect_is(geMegaEnv, "TD")
  expect_equal(dim(geMegaEnvTot), c(45, 4))
  expect_equal(as.numeric(geMegaEnvTot[["megaEnv"]]),
               rep(x = c(2, 1, 1), each = 15))
  expect_equal(levels(geMegaEnvTot[["megaEnv"]]), c("1", "2"))
})

summ <- attr(x = geMegaEnv, which = "sumTab")
test_that("summary is computed correctly", {
  expect_is(summ, "data.frame")
  expect_equal(as.character(summ[["Winning genotype"]]), c("G2", "G2", "G9"))
  expect_equal(summ[["AMMI estimates"]],
               c(111.48629093387, 119.206304548313, 126.941193543844))
})

geMegaEnvMin <- gxeMegaEnv(TD = BLUEs, trait = "t1", method = "min",
                           sumTab = FALSE)
geMegaEnvMinTot <- Reduce(f = rbind, x = geMegaEnvMin)
test_that("option method functions properly", {
  expect_equal(as.numeric(geMegaEnvMinTot[["megaEnv"]]),
               rep(x = c(3, 1, 2), each = 15))
  expect_equal(levels(geMegaEnvMinTot[["megaEnv"]]), c("1", "2", "3"))
  expect_equal(as.character(attr(x = geMegaEnvMin,
                                 which = "sumTab")[["Winning genotype"]]),
               c("G10", "G6", "G7"))
  expect_equal(attr(x = geMegaEnvMin, which = "sumTab")[["AMMI estimates"]],
               c(46.7063274979725, 53.6811189560119, 49.8026925088618))
})

test_that("existing megaEnv in TD is overwritten", {
  expect_warning(gxeMegaEnv(TD = geMegaEnv, trait = "t1", sumTab = FALSE),
                 "TD already contains a column megaEnv")
})

tmp <- capture_output(gxeMegaEnv(TD = BLUEs, trait = "t1", sumTab = FALSE))
test_that("option sumTab succesfully suppresses output", {
  expect_output(gxeMegaEnv(TD = BLUEs, trait = "t1"), "Mega factor")
  expect_equal(tmp, "")
})

test_that("general checks in gxeTable function properly", {
  expect_error(gxeTable(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeTable(geMegaEnv, trait = "t5"),
               "t5 has to be a column in TD")
  expect_error(gxeTable(geMegaEnv, trait = "t1", useYear = TRUE),
               "year has to be a column in TD")
})

test_that("gxeTable functions correctly", {
  ## More random effects than observations, so empty data.frame returned.
  expect_warning(gxeTable(TD = geMegaEnvMin, trait = "t1"),
                 "Empty data.frame returned")
  expect_warning(geTabLm <- gxeTable(TD = geMegaEnv, trait = "t1"),
                 "mega environments that are based on less than 10 trials")
  expect_is(geTabLm, "list")
  expect_length(geTabLm, 2)
  expect_named(geTabLm, c("predictedValue", "standardError"))
  expect_is(geTabLm$predictedValue, "data.frame")
  expect_is(geTabLm$standardError, "data.frame")
  expect_equivalent(geTabLm$predictedValue[1, ],
                    c(79.2416561439773, 79.4864396648177))
  ## This test works fine in RStudio but gives an error when testing on CRAN.
  ## Therefore added a lower tolerance
  expect_equivalent(geTabLm$standardError[1, ],
                    c(6.83020114988906, 6.38884520952826), tolerance = 1e-7)
  skip_on_cran()
  expect_warning(geTabAs <- gxeTable(TD = geMegaEnv, trait = "t1",
                                     engine = "asreml"),
                 "mega environments that are based on less than 10 trials")
  expect_equivalent(geTabAs$predictedValue[1, ],
                    c(79.2091290087036, 79.4263169941458))
  expect_equivalent(geTabAs$standardError[1, ],
                    c(7.20520150684489, 6.96252065379726))
})

## Modify data so it contains a year variable.
geMegaEnvNw2 <- c(geMegaEnv, geMegaEnv)
names(geMegaEnvNw2) <- paste0("E", 1:6)
class(geMegaEnvNw2) <- "TD"
geMegaEnvNw2[["E1"]]$year <- geMegaEnvNw2[["E2"]]$year <-
  geMegaEnvNw2[["E3"]]$year <- 1
geMegaEnvNw2[["E4"]]$year <- geMegaEnvNw2[["E5"]]$year <-
  geMegaEnvNw2[["E6"]]$year <- 2

test_that("option year in gxeTable functions properly", {
  expect_warning(geTab <- gxeTable(TD = geMegaEnvNw2, trait = "t1",
                                   useYear = TRUE),
                 "mega environments that are based on less than 10 trials")
  expect_equivalent(geTab$predictedValue[1, ],
                    c(75.7340703018517, 72.8306094617738))
  expect_equivalent(geTab$standardError[1, ],
                    c(5.69945316015176, 7.83060677272676))
})

test_that("option year in gxeTable functions properly for asreml", {
  skip_on_cran()
  expect_warning(geTab <- gxeTable(TD = geMegaEnvNw2, trait = "t1",
                                   useYear = TRUE, engine = "asreml"),
                 "mega environments that are based on less than 10 trials")
  expect_equivalent(geTab$predictedValue[1, ],
                    c(76.634397082771, 73.7311917123637))
  expect_equivalent(geTab$standardError[1, ],
                    c(5.87848501191572, 8.30017045501152))
})
