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
  expect_s3_class(geVCLm, "varComp")
  expect_named(geVCLm, c("fitMod", "modDat", "trait", "nestingFactor",
                         "useLocYear", "useRegionLocYear", "fullRandVC",
                         "aovFullFixedMod", "engine", "confoundVars",
                         "diagTabs"))
  expect_s4_class(geVCLm$fitMod, "merMod")
  expect_s3_class(geVCLm$modDat, "data.frame")
  expect_null(geVCLm$nestingFactor)
  expect_type(geVCLm$useLocYear, "logical")
  expect_type(geVCLm$useRegionLocYear, "logical")
  expect_s3_class(geVCLm$fullRandVC, "data.frame")
  expect_s3_class(geVCLm$aovFullFixedMod, "anova")
  expect_type(geVCLm$engine, "character")
  expect_type(geVCLm$confoundVars, "character")
  expect_type(geVCLm$diagTabs, "list")
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
  expect_s3_class(geVCAs$fitMod, "asreml")
  expect_s3_class(geVCAs$modDat, "data.frame")
  expect_null(geVCAs$nestingFactor)
  expect_type(geVCAs$useLocYear, "logical")
  expect_type(geVCAs$useRegionLocYear, "logical")
  expect_s3_class(geVCAs$fullRandVC, "data.frame")
  expect_s3_class(geVCAs$aovFullFixedMod, "anova")
  expect_type(geVCAs$engine, "character")
  expect_type(geVCAs$confoundVars, "character")
  expect_type(geVCAs$diagTabs, "list")
})

test_that("option nestingFactor functions correctly", {
  expect_error(gxeVarComp(TD = BLUEs, trait = "t1", engine = "lme4",
                          nestingFactor = "nest"),
               "nest has to be a column in TD")
  ## Produces warning for zero variance component.
  ## Additionaly lme4 gives a related message. Ignoring that here.
  expect_warning(geVCLmNest <- expect_message(gxeVarComp(TD = BLUEs, trait = "t1",
                                                         engine = "lme4",
                                                         nestingFactor = "regime")),
                 "Possible zero variance components")
  expect_equal(geVCLmNest$nestingFactor, "regime")
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
  expect_s3_class(predVCLm, "data.frame")
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
  expect_s3_class(predVCAs, "data.frame")
  expect_named(predVCAs, c("genotype", "predictedValue", "stdError"))
  expect_equal(predVCAs[["predictedValue"]],
               c(74.4223295766225, 68.585100222066, 103.760055840818,
                 81.2827316729689, 75.2942733590616, 96.3375804962255,
                 64.3462150272252, 97.7792057425648, 84.3206060622952,
                 93.7111341197261, 92.5087380574763, 72.4870355565081,
                 68.5207687021667, 70.6451006628306, 101.446118745798))
})

test_that("option predictLevel in predict function functions correctly", {
  predVCLmTr <- predict(geVCLm, predictLevel = "trial")
  expect_s3_class(predVCLmTr, "data.frame")
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
  expect_s3_class(predVCAsTr, "data.frame")
  expect_named(predVCAsTr, c("genotype", "trial", "predictedValue", "stdError"))
  expect_equal(predVCAsTr[["predictedValue"]],
               c(72.6216775379172, 71.4924190545573, 79.1528921373931,
                 66.7844481833606, 65.6551897000008, 73.3156627828366,
                 101.959403802113, 100.830145318753, 108.490618401589,
                 79.4820796342635, 78.3528211509037, 86.0132942337395,
                 73.4936213203563, 72.3643628369964, 80.0248359198322,
                 94.5369284575201, 93.4076699741603, 101.068143056996,
                 62.5455629885198, 61.41630450516, 69.0767775879958,
                 95.9785537038595, 94.8492952204996, 102.509768303335,
                 82.5199540235898, 81.39069554023, 89.0511686230658,
                 91.9104820810207, 90.7812235976609, 98.4416966804967,
                 90.7080860187709, 89.5788275354111, 97.2393006182469,
                 70.6863835178027, 69.5571250344429, 77.2175981172787,
                 66.7201166634613, 65.5908581801015, 73.2513312629373,
                 68.8444486241252, 67.7151901407654, 75.3756632236012,
                 99.6454667070924, 98.5162082237326, 106.176681306568))
})

## VC function

test_that("vc function functions correctly", {
  expect_error(vc(1), "should be an object of class varComp")

  vcVCLm <- vc(geVCLm)
  expect_s3_class(vcVCLm, "data.frame")
  expect_named(vcVCLm, "Component")
  expect_equal(rownames(vcVCLm), c("genotype", "residuals"))
  expect_equal(vcVCLm[["Component"]],
               c(78.7850079727369, 309.69339127849))

  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
  vcVCAs <- vc(geVCAs)
  expect_s3_class(vcVCAs, "data.frame")
  expect_named(vcVCAs, c("Component", "SE"))
  expect_equal(rownames(vcVCAs), c("genotype", "residuals"))
  expect_equal(vcVCAs[["Component"]],
               c(79.9663997429178, 306.149240165981))
})

## herit function

test_that("herit function functions correctly", {
  expect_error(herit(1), "should be an object of class varComp")

  heritVCLm <- herit(geVCLm)
  expect_type(heritVCLm, "double")
  expect_equal(heritVCLm, 0.432846277620034)

  ## Produces warning for zero variance component.
  ## Additionaly lme4 gives a related message. Ignoring that here.
  expect_warning(geVCLmNest <- expect_message(gxeVarComp(TD = BLUEs, trait = "t1",
                                                         engine = "lme4",
                                                         nestingFactor = "regime")),
                 "Possible zero variance components")
  heritVCLmNest <- herit(geVCLmNest)
  expect_equal(heritVCLmNest, 0.604176703866825)

  skip_on_cran()
  skip_on_ci()
  geVCAs <- gxeVarComp(TD = BLUEs, trait = "t1", engine = "asreml")
  heritVCAs <- herit(geVCAs)
  expect_type(heritVCAs, "double")
  expect_equal(heritVCAs, 0.439336846186519)
})

## CRDR function

test_that("CRDR function functions correctly", {
  expect_error(CRDR(geVCLm),
               "CRDR can only be computed when a model is fitted with a nesting")

  ## Produces warning for zero variance component.
  ## Additionaly lme4 gives a related message. Ignoring that here.
  expect_warning(geVCLmNest <- expect_message(gxeVarComp(TD = BLUEs, trait = "t1",
                                                         engine = "lme4",
                                                         nestingFactor = "regime")),
                 "Possible zero variance components")

  geCRDR <- CRDR(geVCLmNest)
  expect_type(geCRDR, "double")
  expect_equal(geCRDR, 1.33852403499764)
})

## correlations function

test_that("correlations function functions correctly", {
  expect_error(correlations(geVCLm),
               "correlations can only be computed when a model is fitted with a nesting")

  ## Produces warning for zero variance component.
  ## Additionaly lme4 gives a related message. Ignoring that here.
  expect_warning(geVCLmNest <- expect_message(gxeVarComp(TD = BLUEs, trait = "t1",
                                                         engine = "lme4",
                                                         nestingFactor = "regime")),
                 "Possible zero variance components")

  geCorr <- correlations(geVCLmNest)
  expect_type(geCorr, "list")
  expect_length(geCorr, 3)
  expect_named(geCorr, c("rScen", "rTrScen", "rTrDiffScen"))
  expect_type(geCorr[[1]], "double")
  expect_type(geCorr[[2]], "double")
  expect_type(geCorr[[3]], "double")
  expect_equivalent(unlist(geCorr), c(1, 0.202803988292456, 0.202803988292456))
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

