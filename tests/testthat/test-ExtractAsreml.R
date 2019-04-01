context("ExtractAsreml")

if (requireNamespace("asreml", quietly = TRUE)) {
  testTD <- createTD(data = testData[testData$field == "E1", ],
                     genotype = "seed", repId = "rep",
                     subBlock = "block", rowId = "Y", colId = "X",
                     rowCoord = "Y", colCoord = "X")

  modelAs <- STRunModel(testTD, design = "res.ibd", traits = "t1",
                        engine = "asreml")

  test_that("the output of extract is of the proper type - asreml", {
    expect_is(STExtract(modelAs, what = "BLUEs"), "list")
    expect_is(STExtract(modelAs), "list")
    expect_length(STExtract(modelAs)[[1]], 23)
    expect_is(STExtract(modelAs, what = c("BLUEs", "BLUPs")), "list")
    expect_length(STExtract(modelAs, what = c("BLUEs", "BLUPs"))[[1]], 2)
  })

  extAs <- STExtract(modelAs)[[1]]
  test_that("BLUEs are computed correctly", {
    expect_is(extAs$BLUEs, "data.frame")
    expect_identical(dim(extAs$BLUEs), c(15L, 2L))
    expect_named(extAs$BLUEs, c("genotype", "t1"))
    expect_equal(extAs$BLUEs$t1,
                 c(91.8559516394161, 73.6623287954585, 74.8843349808429,
                   84.5319987913184, 68.9075264154031, 65.0344131821085,
                   55.3031580266259, 77.0436561256374, 62.1348724815413,
                   67.3983425459108, 108.495613219885, 65.3049181280563,
                   114.94121644052, 78.9340989737003, 76.3592209421872))
  })

  test_that("SE of BLUEs are computed correctly", {
    expect_is(extAs$seBLUEs, "data.frame")
    expect_identical(dim(extAs$seBLUEs), c(15L, 2L))
    expect_named(extAs$seBLUEs, c("genotype", "t1"))
    expect_equal(extAs$seBLUEs$t1, rep(x = 14.5728200113826, times = 15))
  })

  test_that("BLUPs are computed correctly", {
    expect_is(extAs$BLUPs, "data.frame")
    expect_identical(dim(extAs$BLUPs), c(15L, 2L))
    expect_named(extAs$BLUPs, c("genotype", "t1"))
    expect_equal(extAs$BLUPs$t1,
                 c(83.1254685686293, 75.9302650186231, 76.6335759637273,
                   78.3870979895003, 75.3085056773248, 73.8785718154553,
                   68.9076125899214, 75.5226969514441, 72.7178660892812,
                   74.782809512082, 89.2545078476149, 72.9693083410864,
                   90.019070308295, 79.1438169239316, 78.2104770916953))
  })

  test_that("SE of BLUPs are computed correctly", {
    expect_is(extAs$seBLUPs, "data.frame")
    expect_identical(dim(extAs$seBLUPs), c(15L, 2L))
    expect_named(extAs$seBLUPs, c("genotype", "t1"))
    expect_equal(extAs$seBLUPs$t1, rep(x = 8.96361730231537, times = 15))
  })

  test_that("unit errors are computed correctly", {
    expect_is(extAs$ue, "data.frame")
    expect_identical(dim(extAs$ue), c(15L, 2L))
    expect_named(extAs$ue, c("genotype", "t1"))
    expect_equal(extAs$ue$t1, rep(x = 195.688651981688, times = 15))
  })

  test_that("heritability is computed correctly", {
    expect_is(extAs$heritability, "numeric")
    expect_length(extAs$heritability, 1)
    expect_named(extAs$heritability, "t1")
    expect_equivalent(extAs$heritability, 0.345185100576464)
  })

  test_that("varGen is computed correctly", {
    expect_is(extAs$varGen, "numeric")
    expect_length(extAs$varGen, 1)
    expect_named(extAs$varGen, "t1")
    expect_equivalent(extAs$varGen, 99.2921060589769)
  })

  test_that("varErr is computed correctly", {
    expect_is(extAs$varErr, "numeric")
    expect_length(extAs$varErr, 1)
    expect_named(extAs$varErr, "t1")
    expect_equivalent(extAs$varErr, 320.569255027489)
  })

  test_that("fitted values are computed correctly", {
    expect_is(extAs$fitted, "data.frame")
    expect_identical(dim(extAs$fitted), c(30L, 3L))
    expect_named(extAs$fitted, c("genotype", "repId", "t1"))
    expect_equal(extAs$fitted$t1,
                 c(68.7325142281377, 63.8332414120377, 65.1408096443281,
                   95.2835477394968, 48.3415831521691, 61.4693120482354,
                   75.9281618109135, 93.7839958782487, 105.68921935359,
                   124.19321352745, 82.7008157447754, 55.1142370860309,
                   70.963443679784, 88.4283555393336, 72.2011477081544,
                   86.2956532125677, 78.3119310809243, 53.8419769393217,
                   92.7273883030725, 68.5995143159817, 72.7941198083142,
                   79.9243220760605, 67.7916590387072, 71.4567388807611,
                   75.1235098827627, 56.7643391139301, 75.2800017043882,
                   109.956794307189, 107.034432132581, 61.8773220279744))
  })

  test_that("residuals are computed correctly", {
    expect_is(extAs$resid, "data.frame")
    expect_identical(dim(extAs$resid), c(30L, 3L))
    expect_named(extAs$resid, c("genotype", "repId", "t1"))
    expect_equal(extAs$resid$t1,
                 c(-3.73900390681584, 6.66109558076825, -25.5878773515572,
                   9.55404865955619, -10.0401315668725, -2.71192069907498,
                   10.0401315668725, -7.36190314933543, -10.4901112979641,
                   10.490111297964, -8.8681127063347, 8.86811270633465,
                   -6.66109558076828, -9.55404865955656, 20.7840323904322,
                   7.41427077513407, -16.0720951734068, 0.554803295555431,
                   25.5878773515572, 2.71192069907495, 5.81683105357236,
                   -5.81683105357239, -7.41427077513409, 16.0720951734068,
                   -20.7840323904322, -0.554803295555459, 7.3619031493354,
                   4.05951690606135, -4.05951690606138, 3.73900390681578))
  })

  test_that("standardized residuals are computed correctly", {
    expect_is(extAs$stdRes, "data.frame")
    expect_identical(dim(extAs$stdRes), c(30L, 3L))
    expect_named(extAs$stdRes, c("genotype", "repId", "t1"))
    expect_equal(extAs$stdRes$t1,
                 c(-0.322220969016617, 0.574041837408854, -2.20511955611155,
                   0.823351591449873, -0.865241385983749, -0.233708692831061,
                   0.865241385983746, -0.634436236416103, -0.904019870468839,
                   0.904019870468837, -0.764238802846611, 0.764238802846607,
                   -0.574041837408857, -0.823351591449903, 1.79113240419721,
                   0.638949188943159, -1.3850657046552, 0.0478120001911751,
                   2.20511955611155, 0.233708692831059, 0.501284562787259,
                   -0.501284562787262, -0.63894918894316, 1.3850657046552,
                   -1.79113240419721, -0.0478120001911775, 0.6344362364161,
                   0.349842231730744, -0.349842231730747, 0.322220969016612))
  })

  test_that("rMeans are computed correctly", {
    expect_is(extAs$rMeans, "data.frame")
    expect_identical(dim(extAs$rMeans), c(30L, 3L))
    expect_named(extAs$rMeans, c("genotype", "repId", "t1"))
    expect_equal(extAs$rMeans$t1,
                 c(76.292815536098, 68.4662845343126, 62.7335368847438,
                   86.4489757636401, 56.3075860500938, 67.562046837686,
                   83.8948012745056, 92.5987220632557, 85.7264494232054,
                   104.230694382049, 86.485440862549, 58.8982256381372,
                   75.5962544878511, 79.5940275586739, 74.9819671930624,
                   89.7343210251998, 79.9570831587386, 67.9593147643612,
                   90.3207521091556, 74.6920167912245, 71.8939521139256,
                   79.0239220674641, 71.2300760663557, 73.1021349537724,
                   77.9039183264242, 70.8812658977231, 74.0944771044117,
                   91.2281611554149, 88.3062100220531, 69.4378673311319))
  })

  test_that("random effects are computed correctly", {
    expect_is(extAs$ranEf, "data.frame")
    expect_identical(dim(extAs$ranEf), c(15L, 2L))
    expect_named(extAs$ranEf, c("genotype", "t1"))
    expect_equal(extAs$ranEf$t1,
                 c(5.4726918560551, -1.72251169395101, -1.01920074884688,
                   0.734321276926158, -2.34427103524939, -3.77420489711889,
                   -8.74516412265271, -2.13007976113004, -4.93491062329299,
                   -2.8699672004922, 11.6017311350407, -4.68346837148775,
                   12.3662935957208, 1.49104021135746, 0.557700379121089))
  })

  test_that("residuals are computed correctly for genotype random", {
    expect_is(extAs$residR, "data.frame")
    expect_identical(dim(extAs$residR), c(30L, 3L))
    expect_named(extAs$residR, c("genotype", "repId", "t1"))
    expect_equal(extAs$residR$t1,
                 c(-11.2993052147762, 2.02805245849336, -23.1806045919729,
                   18.3886206354129, -18.0061344647972, -8.80465548852551,
                   2.0734921032804, -6.17662933434248, 9.47265863242039,
                   30.452630443365, -12.6527378241084, 5.08412415422836,
                   -11.2939063888354, -0.719720678896905, 18.0032129055242,
                   3.97560296250199, -17.7172472512211, -13.562534529484,
                   27.9945135454742, -3.3805817761678, 6.71699874796091,
                   -4.91643104497604, -10.8526878027826, 14.4266991003954,
                   -23.5644408340937, -14.6717300793484, 8.54742774931188,
                   22.7881500578353, 14.6687052044661, -3.82154139634166))
  })

  test_that("standardized residuals are computed correctly for genotype random", {
    expect_is(extAs$stdResR, "data.frame")
    expect_identical(dim(extAs$stdResR), c(30L, 3L))
    expect_named(extAs$stdResR, c("genotype", "repId", "t1"))
    expect_equal(extAs$stdResR$t1,
                 c(-0.773409428008369, 0.138815162709576, -1.58665491351815,
                   1.25865549227747, -1.23247526218243, -0.602656839137886,
                   0.141925393738444, -0.422774964450603, 0.648379998182383,
                   2.08440705378535, -0.866048534602516, 0.347995693478855,
                   -0.773039891757825, -0.0492630960940438, 1.2322752886935,
                   0.272120166220144, -1.21270164864037, -0.928321863460678,
                   1.91615504644163, -0.231392441229249, 0.459761911095192,
                   -0.336517515905761, -0.742839573439185, 0.987471785848796,
                   -1.61292755266054, -1.00424354886203, 0.585050238123786,
                   1.55979237365422, 1.00403650367579, -0.261575033975524))
  })

  test_that("Waldtest is computed correctly", {
    expect_is(extAs$wald$t1$Ftest, "numeric")
    expect_length(extAs$wald$t1$Ftest, 4)
    expect_equivalent(unlist(extAs$wald$t1$Ftest),
                      c(1.46, 14, 14, 0.244132249447016))
  })

  test_that("CV is computed correctly", {
    expect_is(extAs$CV, "numeric")
    expect_length(extAs$CV, 1)
    expect_named(extAs$CV, "t1")
    expect_equivalent(extAs$CV, 23.4580700451519)
  })

  test_that("rDf is computed correctly", {
    expect_is(extAs$rDf, "integer")
    expect_length(extAs$rDf, 1)
    expect_named(extAs$rDf, "t1")
    expect_equivalent(extAs$rDf, 14)
  })

  test_that("sed is computed correctly", {
    expect_is(extAs$sed$t1, "numeric")
    expect_length(extAs$sed$t1, 3)
    expect_equivalent(extAs$sed$t1,
                      c(18.215842753241, 20.2671887094224, 20.6090797021193))
  })

  test_that("lsd is computed correctly", {
    expect_is(extAs$lsd$t1, "numeric")
    expect_length(extAs$lsd$t1, 3)
    expect_equivalent(extAs$lsd$t1,
                      c(39.0690970463554, 43.4687965454872, 44.2020797953426))
  })

  test_that("correct attributes are added", {
    expect_equal(attr(x = extAs, which = "traits"), "t1")
    expect_equal(attr(x = extAs, which = "design"), "res.ibd")
    expect_equal(attr(x = extAs, which = "engine"), "asreml")
  })

  test_that("calculated values are logically correct", {
    ## Fitted + residuals should match raw data.
    expect_equal(testTD[[1]]$t1 - extAs$fitted$t1 - extAs$resid$t1, rep(0, 30))
    expect_equal(testTD[[1]]$t1 - extAs$rMeans$t1 - extAs$residR$t1, rep(0, 30))
  })
}

