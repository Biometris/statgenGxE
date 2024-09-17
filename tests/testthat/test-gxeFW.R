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
               c(50.0516631185851, 42.1701184913307, 34.5337382287539,
                 16.0842610984795, 7.56166535552369, 2.20832215561466,
                 0.268461475539625, -2.94249396862009, -8.09714897312463,
                 -10.0184273382295, -12.1338406656539, -19.2612066556719,
                 -20.6354998847789, -22.5817246632535, -45.1595111076509))
  expect_equal(est[["SE_Sens"]], rep(x = 28.872707953704, times = 15))
  expect_equal(est[["GenMean"]],
               c(91.7872252465452, 83.4522862345553, 81.6450754872636,
                 101.510790941064, 67.8350451590658, 66.99687610674,
                 60.9764611092266, 77.9806476129257, 94.3724639926346,
                 80.2868026830911, 61.7359582594057, 79.9881044999343,
                 75.6491629852527, 79.1537116210871, 91.2406066244149))
  expect_equal(est[["SE_GenMean"]], rep(x = 9.97742872488619, times = 15))
  expect_equal(est[["MSdeviation"]],
               c(406.747475964132, 968.591074795503, 573.629746805195,
                 1159.71459189021, 722.060188344734, 233.027784089546,
                 288.405924939995, 171.086797099452, 746.29626143719,
                 1072.42047658016, 237.14455249761, 769.684381748083,
                 650.098851153597, 820.783722566764, 139.725726698791))
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
               c(-0.448128709166103, 0.392911552751248, 0.0552171564148552,
                 0.131489446628547, 0.131489446628547, 0.131489446628547,
                 79.2807996719965, 79.9563443298339, 79.685099710811,
                 16.338872847896, 15.1077781545094, 10.1039980734598, 3, 1, 2))
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
