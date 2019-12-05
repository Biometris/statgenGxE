context("gxeFW")

test_that("general checks in gxeFw function properly", {
  expect_error(gxeFw(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeFw(testTD, trait = "t5"),
               "t5 has to be a column in TD")
  expect_error(gxeFw(testTD, trait = "t1", useWt = TRUE),
               "wt has to be a column in TD")
})

geFw0 <- gxeFw(testTD, trait = "t1", maxIter = 30)
test_that("output is of the right class", {
  expect_is(geFw0$estimates, "data.frame")
  expect_is(geFw0$anova, "data.frame")
  expect_is(geFw0$envEffs, "data.frame")
  expect_is(geFw0$TD, "TD")
  expect_is(geFw0$fittedGeno ,"numeric")
  expect_is(geFw0$trait, "character")
  expect_is(geFw0$nGeno, "integer")
  expect_is(geFw0$nEnv, "integer")
  expect_is(geFw0$tol, "numeric")
  expect_is(geFw0$iter, "numeric")
})

est <- geFw0$estimates
test_that("estimates are correct", {
  expect_equal(as.numeric(est[["genotype"]]),
               c(15, 8, 7, 6, 14, 13, 2, 4, 1, 11, 9, 12, 3, 5, 10))
  expect_equal(est[["sens"]],
               c(8.2480292503711, 5.60960266192104, 3.62732300117801,
                 3.09347290553405, 1.57001038415447, 0.870980734809107,
                 0.787095404278508, -0.0536001988633171, -0.231929865723682,
                 -0.295612772451198, -0.534950995967167, -0.767951821629662,
                 -1.37813036252387, -2.44440337065015, -3.09993495443724))
  expect_equal(est[["se_sens"]], rep(x = 2.81089092219227, times = 15))
  expect_equal(est[["genMean"]],
               c(81.4068278048538, 74.8374571665672, 54.2031334320118,
                 97.3247024025305, 71.7251608269665, 94.1532040202909,
                 66.2676045412444, 73.0871656828909, 93.6809432528475,
                 92.2047075941901, 74.4525075789223, 73.2435501712992,
                 81.921276388438, 75.3836554703585, 80.639017569703))
  expect_equal(est[["se_genMean"]], rep(x = 9.73259486220765, times = 15))
  expect_equal(est[["MSdeviation"]],
               c(142.922134663463, 166.778585178951, 48.0878010509886,
                 1897.49132265594, 1037.97082484175, 560.798661088039,
                 422.11953843616, 478.258413915381, 185.829517681951,
                 1076.71319803838, 569.253718630158, 191.840515482879,
                 784.275223041641, 271.634905691247, 691.131887271431))
  expect_equal(est[["rank"]], 1:15)
})

test_that("anova table is correct", {
  expect_equal(as.numeric(as.matrix(geFw0$anova)),
               c(14, 2, 14, 59, 89, 11528.469513222, 1144.51485045288,
                 9254.09383057355, 34100.4249906735, 56027.5031849219,
                 823.462108087286, 572.257425226438, 661.006702183825,
                 577.973304926669, 629.522507695752, 1.42474072949055,
                 0.990110478024662, 1.14366303175142, NA, NA, 0.17069201885813,
                 0.377623515114938, 0.341739661200379, NA, NA))
})

test_that("environmental effects are correct", {
  expect_equal(as.numeric(as.matrix(geFw0$envEffs[, 2:5])),
               c(-0.108308308503234, -4.18529835840507, 4.29360666690831,
                 4.68853685469476, 4.68853685469476, 4.68853685469476,
                 78.8599510406325, 74.7836723680752, 83.2625593719151, 2, 3, 1))
})

test_that("settings for tolerance and maxIter work correctly", {
  expect_error(gxeFw(testTD, trait = "t1", maxIter = 0),
               "a single numerical value greater than or equal to 1")
  expect_error(gxeFw(testTD, trait = "t1", tol = 0),
               "a single numerical value greater than 0")
  expect_warning(gxeFw(testTD, trait = "t1"),
                 "Convergence not achieved in ")
  expect_is(gxeFw(testTD, trait = "t1", maxIter = 30), "FW")
  expect_is(gxeFw(testTD, trait = "t1", tol  = 0.01), "FW")
})

test_that("option sorted sorts estimates correctly", {
  expect_equal(as.numeric(gxeFw(testTD, trait = "t1", maxIter = 30,
                                sorted = "ascending")$estimates[["genotype"]]),
               c(10, 5, 3, 12, 9, 11, 1, 4, 2, 13, 14, 6, 7, 8, 15))
  expect_equal(as.numeric(gxeFw(testTD, trait = "t1", maxIter = 30,
                                sorted = "descending")$estimates[["genotype"]]),
               c(15, 8, 7, 6, 14, 13, 2, 4, 1, 11, 9, 12, 3, 5, 10))
  expect_equal(as.numeric(gxeFw(testTD, trait = "t1", maxIter = 30,
                                sorted = "none")$estimates[["genotype"]]), 1:15)
})

test_that("NA in trait causes no problems", {
  geFw <- gxeFw(testTD, trait = "t4", maxIter = 30)
  expect_is(geFw, "FW")
  expect_equal(dim(geFw$estimates), c(15, 7))
  expect_equal(is.na(c(geFw$TD[["E1"]][["t4"]], geFw$TD[["E2"]][["t4"]],
                       geFw$TD[["E3"]][["t4"]])), is.na(geFw$fittedGeno))
})

test_that("option genotypes functions properly", {
  expect_error(gxeFw(testTD, trait = "t1", genotypes = 1:3),
               "genotypes should be NULL or a character vector")
  expect_error(gxeFw(testTD, trait = "t1", genotypes = "g1"),
               "All genotypes to include should be in TD")
  geFw <- gxeFw(testTD, trait = "t1", genotypes = paste0("G", 1:10))
  expect_is(geFw, "FW")
  expect_equal(dim(geFw$estimates), c(10, 7))
})
