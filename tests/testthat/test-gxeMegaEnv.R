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
               rep(x = c(3, 2, 1), each = 15))
  expect_equal(levels(geMegaEnvTot[["megaEnv"]]), c("1", "2", "3"))
})

summ <- attr(x = geMegaEnv, which = "sumTab")
test_that("summary is computed correctly", {
  expect_is(summ, "data.frame")
  expect_equal(as.character(summ[["Winning genotype"]]), c("G14", "G4", "G5"))
  expect_equal(summ[["AMMI estimates"]],
               c(137.826424971715, 109.450972152511, 127.386419996042))
})

geMegaEnvMin <- gxeMegaEnv(TD = BLUEs, trait = "t1", method = "min",
                           sumTab = FALSE)
geMegaEnvMinTot <- Reduce(f = rbind, x = geMegaEnvMin)
test_that("option method functions properly", {
  expect_equal(as.numeric(geMegaEnvMinTot[["megaEnv"]]),
               rep(x = c(3, 2, 1), each = 15))
  expect_equal(levels(geMegaEnvMinTot[["megaEnv"]]), c("1", "2", "3"))
  expect_equal(as.character(attr(x = geMegaEnvMin,
                                 which = "sumTab")[["Winning genotype"]]),
               c("G12", "G15", "G6"))
  expect_equal(attr(x = geMegaEnvMin, which = "sumTab")[["AMMI estimates"]],
               c(55.0432001091191, 40.7120064386264, 35.3752613055704))
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

## Manipulate geMegaEnv to get proper output.
geMegaEnvNw <- geMegaEnv
geMegaEnvNw[["E3"]][["megaEnv"]] <- 2
test_that("gxeTable functions correctly", {
  ## More random effects than observations, so empty data.frame returned.
  expect_warning(gxeTable(TD = geMegaEnv, trait = "t1"),
                 "Empty data.frame returned")
  expect_warning(geTabLm <- gxeTable(TD = geMegaEnvNw, trait = "t1"),
                 "mega environments that are based on less than 10 trials")
  expect_is(geTabLm, "list")
  expect_length(geTabLm, 2)
  expect_named(geTabLm, c("predictedValue", "standardError"))
  expect_is(geTabLm$predictedValue, "data.frame")
  expect_is(geTabLm$standardError, "data.frame")
  expect_equivalent(geTabLm$predictedValue[1, ],
                    c(80.0426463992245, 80.4692101166694))
  ## This test works fine in RStudio but gives an error when testing on CRAN.
  ## Therefore added a lower tolerance
  expect_equivalent(geTabLm$standardError[1, ],
                    c(5.93814290627681, 9.51137000477223), tolerance = 1e-7)
  skip_on_cran()
  expect_warning(geTabAs <- gxeTable(TD = geMegaEnvNw, trait = "t1",
                                     engine = "asreml"),
                 "mega environments that are based on less than 10 trials")
  expect_equivalent(geTabAs$predictedValue[1, ],
                    c(80.0424619340827, 80.4703103942952))
  expect_equivalent(geTabAs$standardError[1, ],
                    c(6.58334962888207, 9.77290639163001))
})

## Modify data so it contains a year variable.
geMegaEnvNw2 <- c(geMegaEnvNw, geMegaEnvNw)
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
                    c(85.7139211061344, 76.5897665641836))
  expect_equivalent(geTab$standardError[1, ],
                    c(6.72132965673734, 9.57033282833413))
})

test_that("option year in gxeTable functions properly for asreml", {
  skip_on_cran()
  expect_warning(geTab <- gxeTable(TD = geMegaEnvNw2, trait = "t1",
                                   useYear = TRUE, engine = "asreml"),
                 "mega environments that are based on less than 10 trials")
  expect_equivalent(geTab$predictedValue[1, ],
                    c(86.5858674442788, 77.4617649806083))
  expect_equivalent(geTab$standardError[1, ],
                    c(6.98377813814263, 10.2119509172305))
})
