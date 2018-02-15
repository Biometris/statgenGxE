context("gxeFW")

testTD <- createTD(data = testData, genotype = "seed", env = "field",
                   repId = "rep", subBlock = "block", rowId = "Y", colId = "X",
                   rowCoordinates = "Y", colCoordinates = "X")

test_that("settings for tolerance and maxIter work", {
  expect_warning(gxeFw(testTD, trait = "t1"), "Convergence not achieved in ")
  expect_is(gxeFw(testTD, trait = "t1", maxIter = 30), "FW")
  expect_is(gxeFw(testTD, trait = "t1", tol  = 0.01), "FW")
})

geFw <- gxeFw(testTD, trait = "t1", maxIter = 30)
test_that("output is of the right class", {
  expect_is(geFw$estimates, "data.frame")
  expect_is(geFw$anova, "data.frame")
  expect_is(geFw$envEffs, "data.frame")
  expect_is(geFw$TD, "TD")
  expect_is(geFw$fittedGeno ,"numeric")
  expect_is(geFw$trait, "character")
  expect_is(geFw$nGeno, "integer")
  expect_is(geFw$nEnv, "integer")
  expect_is(geFw$tol, "numeric")
  expect_is(geFw$iter, "numeric")
})

est <- geFw$estimates
test_that("estimates are correct", {
  expect_equal(as.numeric(est$genotype),
               c(10, 5, 3, 12, 9, 11, 1, 4, 2, 13, 14, 6, 7, 8, 15))
  expect_equal(est$sens, c(-3.09993495443723, -2.44440337065014, -1.37813036252387,
                           -0.76795182162966, -0.534950995967165, -0.295612772451201,
                           -0.231929865723681, -0.0536001988633181, 0.787095404278507,
                           0.870980734809101, 1.57001038415447, 3.09347290553405,
                           3.62732300117801, 5.60960266192104, 8.24802925037108))
  expect_equal(est$sigmaE, rep(x = 2.81089092219227, times = 15))
  expect_equal(est$genMean, c(80.639017569703, 75.3836554703586, 81.9212763884381,
                              73.2435501712992, 74.4525075789223, 92.2047075941901,
                              93.6809432528475, 73.0871656828909, 66.2676045412444,
                              94.1532040202908, 71.7251608269665, 97.3247024025305,
                              54.2031334320118, 74.8374571665672, 81.4068278048538))
  expect_equal(est$sigma, rep(x = 9.73259486220765, times = 15))
  expect_equal(est$mse, c(691.131887271431, 271.634905691248, 784.275223041641,
                          191.84051548288, 569.253718630158, 1076.71319803839,
                          185.829517681951, 478.258413915381, 422.119538436159,
                          560.798661088037, 1037.97082484175, 1897.49132265594,
                          48.0878010509886, 166.778585178952, 142.922134663462))
})

test_that("anova table is correct", {
  expect_equal(as.numeric(as.matrix(geFw$anova)),
               c(14, 2, 14, 59, 89, 11528.469513222, 1144.51485045288,
                 9254.09383057355, 34100.4249906735, 56027.5031849219,
                 823.462108087286, 572.257425226438, 661.006702183825,
                 577.973304926669, 629.522507695752, 1.42474072949055,
                 0.990110478024662, 1.14366303175142, NA, NA, 0.17069201885813,
                 0.377623515114938, 0.341739661200379, NA, NA))
})

test_that("environmental effects are correct", {
  expect_equal(as.numeric(as.matrix(geFw$envEffs[, 2:5])),
               c(-0.108308308503234, -4.18529835840507, 4.29360666690831,
                 1.37908319436677, 1.37908319436677, 1.37908319436677,
                 78.8599510406325, 74.7836723680752, 83.2625593719151, 2, 3, 1))
})

test_that("option sorted sorts estimates correctly", {
  expect_equal(as.numeric(gxeFw(testTD, trait = "t1", maxIter = 30,
                                sorted = "ascending")$estimates$genotype),
               c(10, 5, 3, 12, 9, 11, 1, 4, 2, 13, 14, 6, 7, 8, 15))
  expect_equal(as.numeric(gxeFw(testTD, trait = "t1", maxIter = 30,
                                sorted = "descending")$estimates$genotype),
               c(15, 8, 7, 6, 14, 13, 2, 4, 1, 11, 9, 12, 3, 5, 10))
  expect_equal(as.numeric(gxeFw(testTD, trait = "t1", maxIter = 30,
                                sorted = "none")$estimates$genotype), 1:15)
})

geFw3 <- gxeFw(testTD, trait = "t3", maxIter = 30)
test_that("NA in trait causes no problems", {
  expect_is(geFw3, "FW")
  expect_identical(dim(geFw3$estimates), c(15L, 6L))
  expect_identical(is.na(geFw3$TD$t3), is.na(geFw3$fittedGeno))
})
