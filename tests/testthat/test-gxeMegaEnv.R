context("gxeMegaEnv")

testTD <- createTD(data = testData, genotype = "seed",
                   trial = "field", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoord = "Y", colCoord = "X")
modelSp <- STRunModel(testTD, design = "rowcol", traits = "t1")
BLUEs <- SSAtoTD(modelSp, what = "BLUEs")

geMegaEnv <- gxeMegaEnv(TD = BLUEs, trait = "t1", sumTab = FALSE)
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

geMegaEnvMin <- gxeMegaEnv(TD = BLUEs, trait = "t1", method = "min",
                           sumTab = FALSE)
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

test_that("combCor helper function funcions correctly", {
  set.seed(1234)
  Xi <- runif(n = 5)
  Yi <- runif(n = 5)
  SXi <- runif(n = 5, max = 0.5)
  SYi <- runif(n = 5, max = 0.5)
  ni <- runif(n = 5, min = 10, max = 20)
  ri <- runif(n = 5)
  r <- combCor(Xi = Xi, Yi = Yi, SXi = SXi, SYi = SYi, ri = ri, ni = ni)
  expect_is(r, "numeric")
  expect_length(r, 1)
  expect_equal(r, 0.146437103999452)
})

test_that("combLocs helper function functions correctly", {
  set.seed(1234)
  Xi <- matrix(runif(n = 9, max = 0.5), nrow = 3,
               dimnames = list(1:3, c("l1", "l2", "l3")))
  SXi <- matrix(runif(n = 9), nrow = 3,
                dimnames = list(1:3, c("l1", "l2", "l3")))
  ammi <- data.frame(genotype = rep(x = paste0("G", 1:10), times = 3),
                     year = rep(1:3, each = 10),
                     l1 = sample(0:1, size = 30, replace = TRUE),
                     l2 = sample(0:1, size = 30, replace = TRUE),
                     l3 = c(sample(0:1, size = 20, replace = TRUE),
                            rep(x = NA, 10)))
  r0 <- by(data = ammi[, c(3:ncol(ammi))], INDICES = ammi$year,
           FUN = cor, use = "pairwise.complete.obs")
  r12 <- combLocs(l1 = "l1", l2 = "l2", ammi = ammi, r0 = r0, Xi = Xi,
                  SXi = SXi)
  expect_is(r12, "numeric")
  expect_length(r12, 1)
  expect_equal(r12, 0.0671135566187331)
  r13 <- combLocs(l1 = "l1", l2 = "l3", ammi = ammi, r0 = r0, Xi = Xi,
                  SXi = SXi)
  expect_equal(r13, -0.0187937478121317)
})
