context("gxeStability")

test_that("general checks in gxeStability function properly", {
  expect_error(gxeStability(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeStability(testTD, trait = c("t1", "t1")),
               "trait has to be a character string of length 1")
  expect_error(gxeStability(testTD, trait = "t5"),
               "t5 has to be a column in TD")
})

geStab <- gxeStability(TD = testTD, trait = "t1")
sup <- geStab$superiority
wri <- geStab$wricke
sta <- geStab$static

test_that("output has the correct structure", {
  expect_is(geStab, "stability")
  expect_is(sup, "data.frame")
  expect_is(sta, "data.frame")
  expect_is(wri, "data.frame")
  expect_identical(geStab$trait, "t1")
})

test_that("superiority is computed correctly", {
  expect_equal(dim(sup), c(15, 3))
  expect_equal(as.numeric(sup[["Genotype"]]),
               c(2, 7, 13, 14, 11, 1, 12, 5, 9, 4, 8, 15, 6, 10, 3))
  expect_equal(sup[["Superiority"]],
               c(2650.29204957547, 2641.97554842648, 2230.66939628004,
                 2228.43156213766, 1787.879596279, 1628.45440542442,
                 1583.88544775107, 1540.50209108773, 1478.13226223008,
                 1383.20034463312, 1249.18177514445, 1195.54492396398,
                 908.409884242796, 864.497039525375, 509.675665244878))
  expect_equal(sup[["Mean"]],
               c(60.9764611092265, 61.7359582594057, 66.99687610674,
                 67.8350451590657, 75.6491629852527, 79.1537116210871,
                 77.9806476129257, 79.9881044999343, 80.2868026830911,
                 81.6450754872636, 83.4522862345553, 91.2406066244148,
                 91.7872252465452, 94.3724639926346, 101.510790941064))
})

test_that("static is computed correctly", {
  expect_equal(dim(sta), c(15, 3))
  expect_equal(as.numeric(sta[["Genotype"]]),
               c(4, 6, 15, 8, 14, 1, 11, 12, 5, 10, 3, 9, 7, 13, 2))
  expect_equal(sta[["Static"]],
               c(503.289300296536, 485.429970860536, 366.941744597168,
                 321.97221439706, 152.542362448676, 116.569131738625,
                 89.4585421269837, 88.7243436702008, 80.4467172774754,
                 66.4221691873904, 46.4243289992708, 42.4810700551934,
                 26.4183911927814, 6.31998874789332, 3.00896819332002))
  expect_equal(sta[["Mean"]],
               c(81.6450754872636, 91.7872252465452, 91.2406066244148,
                 83.4522862345553, 67.8350451590657, 79.1537116210871,
                 75.6491629852527, 77.9806476129257, 79.9881044999343,
                 94.3724639926346, 101.510790941064, 80.2868026830911,
                 61.7359582594057, 66.99687610674, 60.9764611092265))
})

test_that("wricke is computed correctly", {
  expect_equal(dim(wri), c(15, 3))
  expect_equal(as.numeric(wri[["Genotype"]]),
               c(4, 6, 15, 8, 14, 1, 11, 12, 5, 10, 9, 3, 7, 13, 2))
  expect_equal(wri[["Wricke"]],
               c(972.805632845642, 947.433230597357, 759.12178834341,
                 621.541247540926, 311.015799099626, 250.648143023461,
                 188.078360593939, 171.66134793496, 169.17206598567,
                 143.999408512898, 86.9113209375728, 84.1530211193672,
                 60.3173365061717, 13.6368656193455, 7.6276533523883))
  expect_equal(wri[["Mean"]],
               c(81.6450754872636, 91.7872252465452, 91.2406066244148,
                 83.4522862345553, 67.8350451590657, 79.1537116210871,
                 75.6491629852527, 77.9806476129257, 79.9881044999343,
                 94.3724639926346, 80.2868026830911, 101.510790941064,
                 61.7359582594057, 66.99687610674, 60.9764611092265))
})

test_that("option trials functions properly", {
  geStabTr <- gxeStability(TD = testTD, trait = "t1", trials = c("E1", "E2"))
  expect_error(gxeStability(TD = testTD, trait = "t1", trials = "E5"),
               "a character vector defining trials in testTD")
  expect_equal(geStabTr$superiority[["Superiority"]],
               c(2831.58267149612, 2762.09502327719, 2750.20505673928,
                 2421.15454825163, 1891.61485101222, 1773.54921950719,
                 1527.13933756491, 1417.68290034445, 1417.30989947021,
                 1367.35648190921, 1323.70251233194, 1150.64759568343,
                 1097.79517570998, 966.113451036244, 570.09398019407))
})

geStabMin <- gxeStability(TD = testTD, trait = "t1", bestMethod = "min")
test_that("option bestMethod for superiority is computed correctly", {
  supMin <- geStabMin$superiority
  expect_equal(dim(supMin), c(15, 3))
  expect_equal(as.numeric(supMin[["Genotype"]]),
               c(3, 10, 15, 6, 8, 4, 5, 1, 9, 12, 11, 14, 13, 7, 2))
  expect_equal(supMin[["Superiority"]],
               c(1750.83827973778, 1396.13350956294, 1380.45781521972,
                 1329.49473715585, 905.387478755802, 891.42013650905,
                 763.743640669908, 754.368535253886, 749.801240237821,
                 668.697130156506, 614.777532529004, 379.776096309765,
                 311.011548685643, 219.057298956358, 181.898077105332))
  expect_equal(supMin[["Mean"]],
               c(101.510790941064, 94.3724639926346, 91.2406066244148,
                 91.7872252465452, 83.4522862345553, 81.6450754872636,
                 79.9881044999343, 79.1537116210871, 80.2868026830911,
                 77.9806476129257, 75.6491629852527, 67.8350451590657,
                 66.99687610674, 61.7359582594057, 60.9764611092265))
})


test_that("option bestMethod for static is computed correctly", {
  expect_identical(geStabMin$static, geStab$static)
})

test_that("option bestMethod for wricke is computed correctly", {
  expect_identical(geStabMin$wricke, geStab$wricke)
})

test_that("option sort ascending functions properly", {
  geStabAsc <- gxeStability(TD = testTD, trait = "t1", sorted = "ascending")
  expect_identical(geStabAsc$superiority[["genotype"]], rev(sup[["genotype"]]))
  expect_identical(geStabAsc$superiority[["superiority"]],
                   rev(sup[["superiority"]]))
  expect_identical(geStabAsc$static[["genotype"]], rev(sta[["genotype"]]))
  expect_identical(geStabAsc$static[["static"]], rev(sta[["static"]]))
  expect_identical(geStabAsc$wricke[["genotype"]], rev(wri[["genotype"]]))
  expect_identical(geStabAsc$wricke[["wricke"]], rev(wri[["wricke"]]))
})

test_that("option sort none functions properly", {
  geStabNo <- gxeStability(TD = testTD, trait = "t1", sorted = "none")
  expect_equal(as.numeric(geStabNo$superiority[["Genotype"]]), 1:15)
  expect_identical(geStabNo$superiority[["Superiority"]],
                   sup[["Superiority"]][order(sup[["Genotype"]])])
  expect_equal(as.numeric(geStabNo$static[["Genotype"]]), 1:15)
  expect_identical(geStabNo$static$static, sta$static[order(sta[["Genotype"]])])
  expect_equal(as.numeric(geStabNo$wricke[["Genotype"]]), 1:15)
  expect_identical(geStabNo$wricke[["Wricke"]],
                   wri[["Wricke"]][order(wri[["Genotype"]])])
})

test_that("data with missing values is handled correctly", {
  stabMiss <- gxeStability(TD = testTD, trait = "t4")
  expect_is(stabMiss, "stability")
})

