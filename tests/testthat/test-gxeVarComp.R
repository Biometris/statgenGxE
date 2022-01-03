context("gxeVarComp")

test_that("general checks in gxeVarComp function properly", {
  expect_error(gxeVarComp(1, trait = "t1"),
               "TD should be a valid object of class TD")
  expect_error(gxeVarComp(testTD, trait = c("t1", "t1")),
               "trait has to be a character string of length 1")
  expect_error(gxeVarComp(BLUEs, trait = "t5"),
               "t5 has to be a column in TD")
  expect_error(gxeVarComp(BLUEs, trials = "E4", trait = "t1"),
               "a character vector defining trials in BLUEs")
})

geVCLm <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4")
test_that("output structure is of the right class for lme4", {
  expect_is(geVCLm, "varComp")
  expect_named(geVCLm, c("fitMod", "modDat", "trait", "nestingFactor",
                         "useLocYear", "useRegionLocYear", "fullRandVC",
                         "aovFullFixedMod", "engine", "diagTabs"))
  expect_is(geVCLm$fitMod, "merMod")
  expect_is(geVCLm$modDat, "data.frame")
  expect_null(geVCLm$nestingFactor)
  expect_is(geVCLm$useLocYear, "logical")
  expect_is(geVCLm$useRegionLocYear, "logical")
  expect_is(geVCLm$fullRandVC, "data.frame")
  expect_is(geVCLm$aovFullFixedMod, "anova")
  expect_is(geVCLm$engine, "character")
  expect_is(geVCLm$diagTabs, "list")
})

test_that("lme4 model gives the correct output", {
  expect_false(geVCLm$useLocYear)
  expect_named(geVCLm$fullRandVC, c("vcov", "vcovPerc"))
  expect_equal(rownames(geVCLm$fullRandVC), c("trial", "genotype", "residuals"))
  expect_equal(geVCLm$fullRandVC[["vcov"]],
               c(0, 79.9668254460519, 306.14928560952))
  expect_equal(geVCLm$fullRandVC[["vcovPerc"]],
               c(0, 0.207105642982462, 0.792894357017538))
  expect_named(geVCLm$aovFullFixedMod,
               c("Df", "Sum Sq", "Mean Sq", "F value", "Pr(>F)"))
  expect_equal(rownames(geVCLm$aovFullFixedMod),
               c("trial", "genotype", "residuals"))
  expect_equal(geVCLm$aovFullFixedMod[["Sum Sq"]],
               c(513.074183597496, 7644.67782855106, 8671.41494683826))
  expect_equal(geVCLm$aovFullFixedMod[["Mean Sq"]],
               c(256.537091798748, 546.048416325076, 309.693390958509))
  expect_equal(geVCLm$engine, "lme4")
  expect_equal(nrow(geVCLm$diagTabs[[1]]), 0)
})

test_that("output is of the right class for asreml", {
  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
  expect_is(geVCAs$fitMod, "asreml")
  expect_is(geVCAs$modDat, "data.frame")
  expect_null(geVCAs$nestingFactor)
  expect_is(geVCAs$useLocYear, "logical")
  expect_is(geVCAs$fullRandVC, "data.frame")
  expect_is(geVCAs$aovFullFixedMod, "data.frame")
  expect_is(geVCAs$engine, "character")
  expect_is(geVCAs$diagTabs, "list")
})

test_that("option nestingFactor functions correctly", {
  expect_error(gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4",
                          nestingFactor = "nest"),
               "nest has to be a column in TD")
  ## Produces warning for zero variance component.
  ## Ignoring that here.
  # expect_warning(geVCLmNest <- gxeVarComp(TD = BLUEs, trait = "t1",
  #                                         engine = "lme4",
  #                                         nestingFactor = "regime"))
  # expect_equal(geVCLmNest$nestingFactor, "regime")

})

test_that("a warning is printed for likely zero variance components", {
  expect_warning(geVCLm <- gxeVarComp(TD = BLUEs, trait = "t2", engine = "lme4"),
                 "Mean Sum of Squares for genotype smaller")
})

test_that("option diagnostics functions correctly", {
  BLUEs2 <- BLUEs
  ## Add missing value.
  BLUEs2$E1 <- BLUEs2$E1[-1, ]
  expect_output(gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4",
                           diagnostics = TRUE),
                "No missing combinations")
  expect_output(geVCLm2 <- gxeVarComp(TD = BLUEs2, trait = "t1",
                           engine = "lme4", diagnostics = TRUE),
                "1 missing combinations")
  expect_named(geVCLm2$diagTabs[[1]], c("genotype", "trial"))
  expect_equal(as.character(geVCLm2$diagTabs[[1]][["genotype"]]), "G1")
  expect_equal(as.character(geVCLm2$diagTabs[[1]][["trial"]]), "E1")
  ## Add more missing values.
  BLUEs2$E1 <- BLUEs2$E1[-(1:10), ]
  expect_output(gxeVarComp(TD = BLUEs2, trait = "t1",
                           engine = "lme4", diagnostics = TRUE),
                "Printing first 10 missing combinations")
})

## Predict function

test_that("predict function functions correctly", {
  predVCLm <- predict(geVCLm)
  expect_is(predVCLm, "data.frame")
  expect_named(predVCLm, c("genotype", "predictedValue"))
  expect_equal(predVCLm[["predictedValue"]],
               c(79.2972526743991, 76.7659940829781, 92.0192766258885,
                 82.2722003237363, 79.6753627428021, 88.8005912968726,
                 74.9278421249584, 89.425738274789, 83.589545451119,
                 87.6616579859464, 87.1402504519962, 78.4580309618727,
                 76.7380973356425, 77.6592935401848, 91.0158599711691))

  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
  predVCAs <- predict(geVCAs)
  expect_is(predVCAs, "data.frame")
  expect_named(predVCAs, c("genotype", "predictedValue", "stdError"))
  expect_equal(predVCAs[["predictedValue"]],
               c(79.2972526743991, 76.7659940829781, 92.0192766258885,
                 82.2722003237363, 79.6753627428021, 88.8005912968726,
                 74.9278421249584, 89.425738274789, 83.589545451119,
                 87.6616579859464, 87.1402504519962, 78.4580309618727,
                 76.7380973356425, 77.6592935401848, 91.0158599711691))
})

test_that("option predictLevel in predict function functions correctly", {
  predVCLmTr <- predict(geVCLm, predictLevel = "trial")
  expect_is(predVCLmTr, "data.frame")
  expect_named(predVCLmTr, c("genotype", "trial", "predictedValue"))
  expect_equal(predVCLmTr[["predictedValue"]],
               c(77.4966006356938, 74.9653420442727, 90.2186245871831,
                 80.4715482850309, 77.8747107040967, 86.9999392581672,
                 73.127190086253, 87.6250862360836, 81.7888934124136,
                 85.861005947241, 85.3395984132909, 76.6573789231674,
                 74.9374452969371, 75.8586415014795, 89.2152079324637,
                 76.3673421523339, 73.8360835609129, 89.0893661038233,
                 79.3422898016711, 76.7454522207369, 85.8706807748074,
                 71.9979316028932, 86.4958277527238, 80.6596349290537,
                 84.7317474638811, 84.210339929931, 75.5281204398075,
                 73.8081868135772, 74.7293830181196, 88.0859494491039,
                 84.0278152351697, 81.4965566437486, 96.7498391866591,
                 87.0027628845069, 84.4059253035727, 93.5311538576432,
                 79.658404685729, 94.1563008355596, 88.3201080118895,
                 92.3922205467169, 91.8708130127668, 83.1885935226433,
                 81.468659896413, 82.3898561009554, 95.7464225319397))

  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
  predVCAsTr <- predict(geVCAs, predictLevel = "trial")
  expect_is(predVCAsTr, "data.frame")
  expect_named(predVCAsTr, c("genotype", "trial", "predictedValue", "stdError"))
  expect_equal(predVCAsTr[["predictedValue"]],
               c(77.4966006205341, 76.3673421371746, 84.0278152200104,
                 74.9653420188325, 73.8360835354731, 81.4965566183088,
                 90.2186246236927, 89.0893661403333, 96.749839223169,
                 80.4715482819537, 79.3422897985942, 87.00276288143,
                 77.8747106904727, 76.7454522071132, 84.405925289949,
                 86.9999392816045, 85.870680798245, 93.5311538810808,
                 73.1271900533473, 71.9979315699879, 79.6584046528236,
                 87.6250862620599, 86.4958277787004, 94.1563008615362,
                 81.7888934146866, 80.6596349313271, 88.3201080141629,
                 85.8610059660526, 84.7317474826931, 92.3922205655289,
                 85.3395984299848, 84.2103399466253, 91.8708130294611,
                 76.6573789045992, 75.5281204212397, 83.1885935040755,
                 74.9374452713836, 73.8081867880241, 81.4686598708599,
                 75.8586414796673, 74.7293829963079, 82.3898560791436,
                 89.215207964898, 88.0859494815386, 95.7464225643743))
})

## VC function

test_that("vc function functions correctly", {
  expect_error(vc(1), "should be an object of class varComp")

  vcVCLm <- vc(geVCLm)
  expect_is(vcVCLm, "data.frame")
  expect_named(vcVCLm, "Component")
  expect_equal(rownames(vcVCLm), c("genotype", "residuals"))
  expect_equal(vcVCLm[["Component"]],
               c(78.7850079727369, 309.69339127849))

  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
  vcVCAs <- vc(geVCAs)
  expect_is(vcVCAs, "data.frame")
  expect_named(vcVCAs, c("Component", "SE"))
  expect_equal(rownames(vcVCAs), c("genotype", "residuals"))
  expect_equal(vcVCAs[["Component"]],
               c(78.7960647165854, 309.736936464321))
})

## herit function

test_that("herit function functions correctly", {
  expect_error(herit(1), "should be an object of class varComp")

  heritVCLm <- herit(geVCLm)
  expect_is(heritVCLm, "numeric")
  expect_equal(heritVCLm, 0.202804089299665)


  ## Produces warning for zero variance component.
  ## Ignoring that here.
  # expect_warning(geVCLmNest <- gxeVarComp(TD = BLUEs, trait = "t1",
  #                                         engine = "lme4",
  #                                         nestingFactor = "regime"))
  # heritVCLmNest <- herit(geVCLmNest)
  # expect_equal(heritVCLmNest, 0.337218682788646)

  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
  heritVCAs <- herit(geVCAs)
  expect_is(heritVCAs, "numeric")
  expect_equal(heritVCAs, 0.202804046186792)
})

## diagnostics function
test_that("diagnostics function functions correctly", {
  expect_error(diagnostics(1), "should be an object of class varComp")

  BLUEs2 <- BLUEs
  ## Add missing value.
  BLUEs2$E1 <- BLUEs2$E1[-1, ]
  geVCLm <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4")
  expect_output(diagnostics(geVCLm), "No missing combinations")
  geVCLm2 <- gxeVarComp(TD = BLUEs2, trait = "t1", engine = "lme4")
  expect_output(diagnostics(geVCLm2), "1 missing combinations")
  ## Add more missing values.
  BLUEs2$E1 <- BLUEs2$E1[-(1:10), ]
  geVCLm3 <- gxeVarComp(TD = BLUEs2, trait = "t1", engine = "lme4")
  expect_output(diagnostics(geVCLm3), "11 missing combinations")
})


