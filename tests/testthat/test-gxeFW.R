context("gxeFW")

test_that("general checks in gxeFw function properly", {
  expect_error(gxeFw(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeFw(testTD, trait = c("t1", "t1")),
               "trait has to be a character string of length 1")
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
  expect_is(geFw0$fittedGeno ,"data.frame")
  expect_is(geFw0$trait, "character")
  expect_is(geFw0$nGeno, "integer")
  expect_is(geFw0$nEnv, "integer")
  expect_is(geFw0$tol, "numeric")
  expect_is(geFw0$iter, "numeric")
})

est <- geFw0$estimates
test_that("estimates are correct", {
  expect_equal(as.numeric(est[["Genotype"]]),
               c(6, 8, 4, 3, 14, 13, 2, 12, 10, 9, 7, 5, 11, 1, 15))
  expect_equal(est[["Sens"]],
               c(62.3138747318992, 52.5013357473576, 42.9936348166129,
                 20.0247179128608, 9.41449781293969, 2.74939616289115,
                 0.334277144158371, -3.66361671330498, -10.0806564687019,
                 -12.4729518166929, -15.1064789550751, -23.9800680568997,
                 -25.6910435690698, -28.1138405697322, -56.2230781792431))
  expect_equal(est[["SE_Sens"]], rep(x = 35.9459727446209, times = 15))
  expect_equal(est[["GenMean"]],
               c(91.7872252465452, 83.4522862345553, 81.6450754872636,
                 101.510790941064, 67.8350451590658, 66.99687610674,
                 60.9764611092266, 77.9806476129257, 94.3724639926346,
                 80.2868026830911, 61.7359582594057, 79.9881044999343,
                 75.6491629852527, 79.1537116210871, 91.2406066244149))
  expect_equal(est[["SE_GenMean"]], rep(x = 9.97742872488619, times = 15))
  expect_equal(est[["MSdeviation"]],
               c(406.745164521851, 968.590479155985, 573.634227710133,
                 1159.71455627792, 722.05950064614, 233.027744798405,
                 288.405921397104, 171.086587643299, 746.296717901325,
                 1072.4200985057, 237.144572332745, 769.683832455474,
                 650.098279945242, 820.784587275549, 139.725285849686))
  expect_equal(est[["Rank"]], 1:15)
})

test_that("anova table is correct", {
  expect_equal(as.numeric(as.matrix(geFw0$anova)),
               c(2, 14, 14, 59, 89, 9.55053113096801, 11901.6249899836,
                 6771.97406517297, 35837.6702256662, 54520.8198119537,
                 4.77526556548401, 850.116070713112, 483.71243322664,
                 607.418139418072, 612.593481033188, 0.00786157879654184,
                 1.39955660778842, 0.79634176498261, NA, NA, 0.992170281738418,
                 0.182343051582967, 0.669153553775784, NA, NA))
})

test_that("environmental effects are correct", {
  expect_equal(as.numeric(as.matrix(geFw0$envEffs[, 2:6])),
               c(-0.359948234376616, 0.315591707750728, 0.0443565266258885,
                 0.491548965248676, 0.491548965248676, 0.491548965248676,
                 79.2808016613818, 79.9563378658385, 79.6851041854211,
                 16.3388729087687, 15.107650810596, 10.1040248456418, 3, 1, 2))
})

test_that("settings for tolerance and maxIter work correctly", {
  expect_error(gxeFw(testTD, trait = "t1", maxIter = 0),
               "a single numerical value greater than or equal to 1")
  expect_error(gxeFw(testTD, trait = "t1", tol = 0),
               "a single numerical value greater than 0")
  expect_warning(gxeFw(testTD, trait = "t1", maxIter = 2),
                 "Convergence not achieved in ")
  expect_is(gxeFw(testTD, trait = "t1", maxIter = 30), "FW")
  expect_is(gxeFw(testTD, trait = "t1", tol  = 0.01), "FW")
})

test_that("option sorted sorts estimates correctly", {
  expect_equal(as.numeric(gxeFw(testTD, trait = "t1", maxIter = 30,
                                sorted = "ascending")$estimates[["Genotype"]]),
               c(15, 1, 11, 5, 7, 9, 10, 12, 2, 13, 14, 3, 4, 8, 6))
  expect_equal(as.numeric(gxeFw(testTD, trait = "t1", maxIter = 30,
                                sorted = "descending")$estimates[["Genotype"]]),
               c(6, 8, 4, 3, 14, 13, 2, 12, 10, 9, 7, 5, 11, 1, 15))
  expect_equal(as.numeric(gxeFw(testTD, trait = "t1", maxIter = 30,
                                sorted = "none")$estimates[["Genotype"]]), 1:15)
})

test_that("NA in trait causes no problems", {
  geFw <- gxeFw(testTD, trait = "t4", maxIter = 50)
  expect_is(geFw, "FW")
  expect_equal(dim(geFw$estimates), c(15, 7))
  expect_equal(sum(is.na(geFw$fittedGeno[["fittedValue"]])), 0)
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

test_that("a user friendly warning is given for genotypes in single trial", {
  testTD2 <- testTD
  testTD2$E2 <- testTD2$E2[testTD2$E2[["genotype"]] != "G1", ]
  testTD2$E3 <- testTD2$E3[testTD2$E3[["genotype"]] != "G1", ]
  expect_warning(gxeFw(testTD2, trait = "t1"),
                 c("following genotypes are present in only one trial"))
  expect_silent(gxeFw(testTD2, trait = "t1", genotypes = paste0("G", 2:15)))
})
