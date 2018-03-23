context("ExtractAsreml")

testTD <- createTD(data = testData[testData$field == "E1", ],
                   genotype = "seed", repId = "rep",
                   subBlock = "block", rowId = "Y", colId = "X",
                   rowCoordinates = "Y", colCoordinates = "X")

modelAs <- STRunModel(testTD, design = "res.ibd", traits = "t1",
                      engine = "asreml")

test_that("the output of extract is of the proper type - asreml", {
  expect_is(STExtract(modelAs, what = "BLUEs"), "list")
  expect_is(STExtract(modelAs), "list")
  expect_length(STExtract(modelAs)[[1]], 18)
  expect_is(STExtract(modelAs, what = c("BLUEs", "BLUPs")), "list")
  expect_length(STExtract(modelAs, what = c("BLUEs", "BLUPs"))[[1]], 2)
})

extAs <- STExtract(modelAs)[[1]]
test_that("BLUEs are computed correctly", {
  expect_is(extAs$BLUEs, "data.frame")
  expect_identical(dim(extAs$BLUEs), c(15L, 2L))
  expect_equal(colnames(extAs$BLUEs), c("genotype", "t1"))
  expect_equal(extAs$BLUEs$t1, c(91.8559516394161, 73.6623287954585, 74.8843349808429,
                                 84.5319987913184, 68.9075264154031, 65.0344131821085,
                                 55.3031580266259, 77.0436561256374, 62.1348724815413,
                                 67.3983425459108, 108.495613219885, 65.3049181280563,
                                 114.94121644052, 78.9340989737003, 76.3592209421872))
})

test_that("SE of BLUEs are computed correctly", {
  expect_is(extAs$seBLUEs, "data.frame")
  expect_identical(dim(extAs$seBLUEs), c(15L, 2L))
  expect_equal(colnames(extAs$seBLUEs), c("genotype", "t1"))
  expect_equal(extAs$seBLUEs$t1, rep(x = 14.5728200113826, times = 15))
})

test_that("BLUPs are computed correctly", {
  expect_is(extAs$BLUPs, "data.frame")
  expect_identical(dim(extAs$BLUPs), c(15L, 2L))
  expect_equal(colnames(extAs$BLUPs), c("genotype", "t1"))
  expect_equal(extAs$BLUPs$t1, c(83.1254685686293, 75.9302650186231, 76.6335759637273,
                                 78.3870979895003, 75.3085056773248, 73.8785718154553,
                                 68.9076125899214, 75.5226969514441, 72.7178660892812,
                                 74.782809512082, 89.2545078476149, 72.9693083410864,
                                 90.019070308295, 79.1438169239316, 78.2104770916953))
})

test_that("SE of BLUPs are computed correctly", {
  expect_is(extAs$seBLUPs, "data.frame")
  expect_identical(dim(extAs$seBLUPs), c(15L, 2L))
  expect_equal(colnames(extAs$seBLUPs), c("genotype", "t1"))
  expect_equal(extAs$seBLUPs$t1, rep(x = 8.96361730231537, times = 15))
})

test_that("unit errors are computed correctly", {
  expect_is(extAs$ue, "data.frame")
  expect_identical(dim(extAs$ue), c(15L, 2L))
  expect_equal(colnames(extAs$ue), c("genotype", "t1"))
  expect_equal(extAs$ue$t1, rep(x = 195.688651981688, times = 15))
})

test_that("heritability is computed correctly", {
  expect_is(extAs$heritability, "numeric")
  expect_length(extAs$heritability, 1)
  expect_equal(names(extAs$heritability), "t1")
  expect_equal(unname(extAs$heritability), 0.345185100576464)
})

test_that("varGen is computed correctly", {
  expect_is(extAs$varGen, "numeric")
  expect_length(extAs$varGen, 1)
  expect_equal(names(extAs$varGen), "t1")
  expect_equal(unname(extAs$varGen), 99.2921060589769)
})

test_that("varErr is computed correctly", {
  expect_is(extAs$varErr, "numeric")
  expect_length(extAs$varErr, 1)
  expect_equal(names(extAs$varErr), "t1")
  expect_equal(unname(extAs$varErr), 320.569255027489)
})

test_that("fitted values are computed correctly", {
  expect_is(extAs$fitted, "data.frame")
  expect_identical(dim(extAs$fitted), c(30L, 3L))
  expect_equal(colnames(extAs$fitted), c("genotype", "repId", "t1"))
  expect_equal(extAs$fitted$t1, c(95.2835477394974, 86.2956532125678, 75.9281618109139,
                                  70.9634436797838, 109.956794307189, 68.7325142281375,
                                  124.193213527451, 92.7273883030729, 79.9243220760603,
                                  75.1235098827624, 78.3119310809241, 93.7839958782488,
                                  82.7008157447758, 68.5995143159815, 56.7643391139298,
                                  88.4283555393344, 67.791659038707, 48.3415831521687,
                                  63.8332414120378, 107.034432132581, 61.8773220279746,
                                  105.68921935359, 65.1408096443277, 72.7941198083143,
                                  72.2011477081546, 71.4567388807612, 75.280001704388,
                                  55.1142370860305, 61.4693120482355, 53.841976939322))
})

test_that("residuals are computed correctly", {
  expect_is(extAs$resid, "data.frame")
  expect_identical(dim(extAs$resid), c(30L, 3L))
  expect_equal(colnames(extAs$resid), c("genotype", "repId", "t1"))
  expect_equal(extAs$resid$t1, c(9.55404865955563, 7.41427077513394, 10.0401315668721,
                                 -6.66109558076809, 4.05951690606165, -3.73900390681565,
                                 10.4901112979639, 25.5878773515568, -5.8168310535722,
                                 -20.7840323904319, -16.0720951734066, -7.36190314933556,
                                 -8.8681127063351, 2.71192069907514, -0.554803295555175,
                                 -9.55404865955744, -7.41427077513389, -10.0401315668721,
                                 6.66109558076813, -4.05951690606162, 3.73900390681564,
                                 -10.4901112979639, -25.5878773515568, 5.81683105357222,
                                 20.7840323904319, 16.0720951734067, 7.3619031493356,
                                 8.86811270633509, -2.71192069907509, 0.554803295555182))
})

test_that("standardized residuals are computed correctly", {
  expect_is(extAs$stdRes, "data.frame")
  expect_identical(dim(extAs$stdRes), c(30L, 3L))
  expect_equal(colnames(extAs$stdRes), c("genotype", "repId", "t1"))
  expect_equal(extAs$stdRes$t1, c(0.823351591449851, 0.638949188943169, 0.86524138598374,
                                  -0.574041837408859, 0.349842231730781, -0.322220969016611,
                                  0.904019870468856, 2.20511955611159, -0.501284562787263,
                                  -1.79113240419724, -1.38506570465523, -0.634436236416135,
                                  -0.76423880284667, 0.233708692831082, -0.0478120001911546,
                                  -0.823351591450006, -0.638949188943164, -0.86524138598374,
                                  0.574041837408862, -0.349842231730779, 0.32222096901661,
                                  -0.904019870468852, -2.20511955611159, 0.501284562787264,
                                  1.79113240419724, 1.38506570465523, 0.634436236416138,
                                  0.76423880284667, -0.233708692831079, 0.0478120001911552))
})

test_that("rMeans are computed correctly", {
  expect_is(extAs$rMeans, "data.frame")
  expect_identical(dim(extAs$rMeans), c(30L, 3L))
  expect_equal(colnames(extAs$rMeans), c("genotype", "repId", "t1"))
  expect_equal(extAs$rMeans$t1, c(86.4489757636404, 89.7343210252, 83.8948012745057,
                                  75.5962544878507, 91.2281611554156, 76.2928155360976,
                                  104.230694382051, 90.3207521091562, 79.023922067464,
                                  77.9039183264238, 79.9570831587384, 92.5987220632562,
                                  86.4854408625493, 74.692016791224, 70.8812658977221,
                                  79.5940275586746, 71.2300760663556, 56.3075860500929,
                                  68.4662845343125, 88.3062100220543, 69.4378673311317,
                                  85.7264494232064, 62.7335368847434, 71.8939521139257,
                                  74.9819671930626, 73.1021349537726, 74.0944771044117,
                                  58.8982256381365, 67.5620468376858, 67.9593147643609))
})

test_that("random effects are computed correctly", {
  expect_is(extAs$ranEf, "data.frame")
  expect_identical(dim(extAs$ranEf), c(15L, 2L))
  expect_equal(colnames(extAs$ranEf), c("genotype", "t1"))
  expect_equal(extAs$ranEf$t1, c(5.4726918560551, -1.72251169395101, -1.01920074884688,
                                 0.734321276926158, -2.34427103524939, -3.77420489711889,
                                 -8.74516412265271, -2.13007976113004, -4.93491062329299,
                                 -2.8699672004922, 11.6017311350407, -4.68346837148775,
                                 12.3662935957208, 1.49104021135746, 0.557700379121089))
})

test_that("Waldtest is computed correctly", {
  expect_is(extAs$wald$t1$Ftest, "numeric")
  expect_length(extAs$wald$t1$Ftest, 4)
  expect_equal(unname(unlist(extAs$wald$t1$Ftest)), c(1.46, 14, 14, 0.244132249447016))
})

test_that("CV is computed correctly", {
  expect_is(extAs$CV, "numeric")
  expect_length(extAs$CV, 1)
  expect_equal(names(extAs$CV), "t1")
  expect_equal(unname(extAs$CV), 23.4580700451519)
})

test_that("rDf is computed correctly", {
  expect_is(extAs$rDf, "integer")
  expect_length(extAs$rDf, 1)
  expect_equal(names(extAs$rDf), "t1")
  expect_equal(unname(extAs$rDf), 14)
})

test_that("sed is computed correctly", {
  expect_is(extAs$sed$t1, "numeric")
  expect_length(extAs$sed$t1, 3)
  expect_equal(unname(extAs$sed$t1),
               c(18.215842753241, 20.2671887094224, 20.6090797021193))
})

test_that("lsd is computed correctly", {
  expect_is(extAs$lsd$t1, "numeric")
  expect_length(extAs$lsd$t1, 3)
  expect_equal(unname(extAs$lsd$t1),
               c(39.0690970463554, 43.4687965454872, 44.2020797953426))
})

test_that("correct attributes are added", {
  expect_equal(attr(x = extAs, which = "traits"), "t1")
  expect_equal(attr(x = extAs, which = "design"), "res.ibd")
  expect_equal(attr(x = extAs, which = "engine"), "asreml")
})

